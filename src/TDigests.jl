module TDigests

export Cluster, TDigest, mergedata!, insertdata, mergebuffer

include("scales.jl")

"""
    Cluster(x [, count::Int=1])
    
Initialize a `Cluster` of values included in a `TDigest`.

A `Cluster` is characterized by its centroid and the count of values included in
it. Those parameters can be queried by the functions `centroid` and `count`.

If only a value `x` is passed to the constructor of `Cluster`, a singleton
cluster is created, with a centroid equal to `Float64(x)`.
"""
struct Cluster # tbd: support different AbstractFloat types
    centroid::Float64
    count::Int
end
Cluster(x::T) where {T<:Real} = Cluster(Float64(x), 1)

centroid(c::Cluster)  = c.centroid
count(c::Cluster) = c.count

"""
    union(c1::Cluster, c2::Cluster)
    c1 ∪ c2
    
Construct a union of two `Cluster`s of a `TDigest`.

The centroid of the new `Cluster` is the weighted average of `centroid(c1)` and
`centroid(c2)`, and its count is `count(c1) + count(c2)`.
"""
function Base.union(c1::Cluster, c2::Cluster)
    count = c1.count + c2.count
    Cluster((c1.centroid*c1.count + c2.centroid*c2.count)/count, count)
end

# Overloading `Base` functions for `Cluster`
for comp in (:isless, :isequal)
    @eval Base.$comp(c1::Cluster, c2::Cluster) = $comp(c1.centroid, c2.centroid)
end

Base.zero(T::Type{Cluster}) = Cluster(zero(Float64),0)



"""
    TDigest(δ::Int [, k::Function=k1, kinv::Function=k1_inv, buffersize::Int=10δ])
    
Construct an empty `TDigest` with compression parameter `δ`.

A `TDigest` is a data structure composed by a maximum of `⌈δ⌉` clusters that are
arranged in a sorted vector (see [`clusters`](@ref)). It also contains a buffer
that is used to collect new data to be included in the digest through a
progressive merging algorithm (see [`buffer](@ref)). By default the amount
of space allocated for the buffer is 10 times the compression parameter `δ`.

The merging algorithm uses the scale function `k` and its inverse `kinv` to
define what values of the buffer have to be merged with previous centroids.
Cf. the functions [`TDigests.k0`](@ref) and [`TDigests.k1`](@ref). 

## References:

[1] Dunning, T. & Ertl, O. "Computing extremely accurate quantiles using t-digests",
arXiv: 1902.04023v1, 2019.
"""
struct TDigest{F<:Function, FI<:Function}
    clusters::Vector{Cluster}
    buffer::Vector{Cluster}
    δ::Int
    k::F
    kinv::FI
    # actual clusters: .clusters[_limits[1]:_limits[2]]
    # actual buffer: .buffer[1:_limits[3]]
    _limits::Vector{Int}
    # extrema for more accurate estimation of extreme quantiles
    _extrema::Vector{Float64}
    # Inner constructors sets empty centroid and buffer vectors
    function TDigest(δ::Int, k::Function=k1, kinv::Function=k1_inv, buffersize::Int=10δ)
        new{typeof(k),typeof(kinv)}(zeros(Cluster, δ),
                                    zeros(Cluster, δ+buffersize),
                                    δ, k, kinv,
                                    [1, 0, 0], zeros(Float64,2))
    end
end
function TDigest(clusters, buffer, δ, k, kinv)
    nc = length(clusters)
    if nc > δ
        throw(ArgumentError("`length(clusters)` must be greater or equal than `δ`"))
    end
    td = TDigest(δ, k, kinv, length(clusters)+length(buffersize))
    td.clusters .= clusters
    td.buffer   .= buffer
    td._limits  .= [1, nc, nc+length(buffer)]
    td._extrema .= [extrema(centroid.(clusters))...]
end
# Other constructors tbd

"""
    clusters(td::TDigest)
    
Extract a view as a sorted vector of the clusters contained in `td`.
"""
clusters(td::TDigest) = view(td.clusters, td._limits[1]:td._limits[2])

"""
    buffer(td::TDigest)
    
Extract a view of the clusters buffered to be merged in `td`.
"""
buffer(td::TDigest)    = view(td.buffer, 1:td._limits[3])

# Overloading methods for `AbstractVector`.
Base.length(td::TDigest) = td._limits[2] - td._limits[1] + 1
for f in (:getindex, :view, :iterate)
    @eval Base.$f(td::TDigest, args...) = $f(clusters(td), args...)
end

## Functions defined in Dunning's and Ertl's paper, but not really used here

leftweight(td::TDigest, i) = sum(count, view(clusters(td), 1:i-1))
rightweight(td::TDigest, i) = sum(count, view(clusters(td), 1:i))

function ksize(td::TDigest, i)
    n = length(td)
    qleft = leftweight(td, i)/n
    qright = qleft + td[i].count/n
    return td.k(qright) - td.k(qleft)
end


## Merging algorithm ##
const MergeForward = false
const MergeReverse = true

"""
    mergedata!(td::TDigest, newdata)
    
Merge progressively the items of the iterator `newdata` into the t-digest `td`.
"""
function mergedata!(td::TDigest, newdata)
    state = ()                  # for the initial iteration 
    finished = isempty(newdata) # ensure that there are data to iterate on
    forward = true              # start with forward-merging
    while !finished
        n = length(td)
        td.buffer[1:n] .= clusters(td) # fill the buffer with existing clusters
        # update `state` with the last successful iteration
        # until the buffer is full or there is nothing else to iterate (`finished==true`)
        state, finished = insertbuffer!(td, n+1, newdata, state)
        # sort and merge, and repeat if not finished, with reversed merging direction
        b = buffer(td)
        sort!(b)
        if length(td) == 0
            td._extrema[1] = centroid(b[1])
            td._extrema[2] = centroid(b[end])
        else
            (centroid(b[1]) < td._extrema[1]) && (td._extrema[1] = centroid(b[1]))
            (centroid(b[end]) > td._extrema[2]) && (td._extrema[2] = centroid(b[end]))
        end
        if forward
            mergebuffer!(td, Val(MergeForward))
        else
            mergebuffer!(td, Val(MergeReverse))
        end
        forward = !forward
    end
    td._limits[3] = 0
    return nothing
end

"""
    insertbuffer!(td::TDigest, n, newdata, state) -> state, finished::Bool
    
Insert into the buffer of `td` as many items from `newdata` as possible, starting
at the index `n` of the buffer and at the position of the iterator marked by `state`.

The function stops when the buffer is full or when the end of `newdata` is reached.
It returns the state after the last successful iteration, and a `Bool` indicating if 
`newdata` is finished before filling the buffer, in order to set up subsequent
reading operations of `newdata` if necessary.
"""
function insertbuffer!(td::TDigest, n, newdata, state)
    iteration = iterate(newdata, state...)
    # Recursive insertion of iterations into the buffer (function barrier)
    # return the last position inserted, and the state after last successful iteration
    n, state = _bufferiteration!(td, n, iteration, newdata, state)
    td._limits[3] = n
    finished = (n < length(td.buffer))
    return (state, finished)
end

# methods of _bufferiteration (single insertion of data included in `iteration`)

function _bufferiteration!(td, n, iteration, newdata, _) # for successful iterations
    value, state = iteration
    td.buffer[n] = Cluster(value)
    if n < length(td.buffer) # condition for recursion
        iteration = iterate(newdata, state)
        n, state = _bufferiteration!(td, n+1, iteration, newdata, state)
    end
    return (n, state)
end

# for failed iteration (the iteration reached its end)
function _bufferiteration!(td, n, iteration::Nothing, newdata, previous_state)
    return (n-1, previous_state)
end
    
_qlimit(td::TDigest, q0, order::Val{MergeForward}) = td.kinv(td.k(q0, td.δ) + 1.0, td.δ)
_qlimit(td::TDigest, q0, order::Val{MergeReverse}) = td.kinv(td.k(q0, td.δ) - 1.0, td.δ)

function mergebuffer!(td::TDigest, mergeorder::Val{MergeForward})
    bf = buffer(td)
    s = sum(count, bf)
    q0 = 0.0
    qlimit = _qlimit(td, q0, mergeorder)
    σ = bf[1]
    n = 1
    for x ∈ view(bf, 2:length(bf))
        q = q0 + (σ.count + x.count)/s
        if q ≤ qlimit
            σ = σ ∪ x
        else
            td.clusters[n] = σ
            n += 1
            q0 += σ.count/s
            qlimit = _qlimit(td, q0, mergeorder)
            σ = x
        end
    end
    td.clusters[n] = σ
    td._limits .= [1, n, 0]
    return n #length(td)
end

function mergebuffer!(td::TDigest, mergeorder::Val{MergeReverse})
    bf = buffer(td)
    s = sum(count, bf)
    q0 = 1.0
    qlimit = _qlimit(td, q0, mergeorder)
    σ = bf[end]
    n = td.δ
    for x ∈ view(bf, length(bf)-1:-1:1)
        q = q0 - (σ.count + x.count)/s
        if q ≥ qlimit
            σ = σ ∪ x
        else
            td.clusters[n] = σ
            n -= 1
            q0 -= σ.count/s
            qlimit = _qlimit(td, q0, mergeorder)
            σ = x
        end
    end
    td.clusters[n] = σ
    td._limits .= [n, td.δ, 0]
    return n #length(td)
end



function quantile(td::TDigest, φ::Real)
    (φ ≈ 0.0) && return td._extrema[1]
    (φ ≈ 1.0) && return td._extrema[2]
    nc = length(td)
    cc = clusters(td)
    ranks = cumsum(count.(cc))
    φn = φ*ranks[end]
    k = searchsortedfirst(ranks, φn)
    if k == 1      # hit first cluster
        return interpolate_centroid_first(φn, td._extrema[1], cc[1], cc[2])
    elseif k ≥ nc  # hit last cluster
        return interpolate_centroid_first(φn, cc[nc-1], cc[nc], td._extrema[2], ranks[end])
    else
        return interpolate_centroid(φn, ranks[k], cc[k-1], cc[k], cc[k+1])
    end
end

function interpolate_centroid(φn, r, c₋, c, c₊)
    w₋, w, w₊ = count.((c₋, c, c₊))
    # estimate ranks in the midpoints of the lower and upper clusters
    r_lower = (w₋ == 1) ? float(r - w) : float(r - w - w₋/2)
    r_upper = (w₊ == 1) ? float(r + 1) : float(r + w₊/2)
    x_lower, x, x_upper = centroid.((c₋, c, c₊))
    # fix upper or lower values for the interpolation
    if w == 1              # for singleton middle cluster
        r_upper = float(r)
        x_upper = x
    elseif r - φn > (w/2)  # for target in the lower half
        r_upper = float(r - w/2)
        x_upper = x
    else                   # for target in the upper half
        r_lower = float(r - w/2)
        x_lower = x
    end
    _interpolate(φn, (r_lower, r_upper), (x_lower, x_upper))
end

function interpolate_centroid_first(φn, minvalue, c, c₊)
    (φn ≤ 1) && return minvalue
    w, w₊ = count.((c, c₊))
    if w == 1      # singleton first cluster
        return centroid(c)
    elseif w == 2  # twin first cluster
        x_lower = minvalue
        x_upper = 2centroid(c) - minvalue
        return _interpolate(φn, (1, 2), (x_lower, x_upper))
    else
        # add singleton with minimum value and add 1 to the ranks
        return interpolate_centroid(φn+1, w+1, Cluster(minvalue), c, c₊)
    end
end

function interpolate_centroid_last(φn, c₋, c, maxvalue, n)
    φn = n - φn + 1
    interpolate_centroid_first(φn, maxvalue, c, c₋)
end

function _interpolate(x, extrema_x, extrema_y)
    t = (x - extrema_x[1])/(extrema_x[2] - extrema_x[1])
    return extrema_y[1] + (extrema_y[2] - extrema_y[1])*t
end

end # module
