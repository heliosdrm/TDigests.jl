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
# struct Cluster # tbd: support different AbstractFloat types
#    centroid::Float64
#    count::Int
# end
# Cluster(x::T) where {T<:Real} = Cluster(Float64(x), 1)

const Cluster = NamedTuple{(:centroid, :count),Tuple{Float64,Int64}}

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
    # Cluster((c1.centroid*c1.count + c2.centroid*c2.count)/count, count)
    (centroid=(c1.centroid*c1.count + c2.centroid*c2.count)/count, count=count)
end

# Overloading `Base` functions for `Cluster`
for comp in (:isless, :isequal)
    @eval Base.$comp(c1::Cluster, c2::Cluster) = $comp(c1.centroid, c2.centroid)
end

# Base.zero(T::Type{Cluster}) = Cluster(zero(Float64),0)
Base.zero(T::Type{Cluster}) = (centroid=0.0, count=0)


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
mutable struct TDigest{F<:Function, FI<:Function}
    clusters::Vector{Cluster}
    buffer::Vector{Cluster}
    δ::Int
    k::F
    kinv::FI
    # actual clusters: .clusters[_limits[1]:_limits[2]]
    # actual buffer: .buffer[1:_limits[3]]
    _c_first::Int
    _c_last::Int
    _b_first::Int
    _b_last::Int
    # extrema for more accurate estimation of extreme quantiles
    _min::Float64
    _max::Float64
    # Inner constructors sets empty centroid and buffer vectors
    function TDigest(δ::Int, k::Function=k1, kinv::Function=k1_inv, buffersize::Int=10δ)
        new{typeof(k),typeof(kinv)}(zeros(Cluster, δ),
                                    zeros(Cluster, δ+buffersize),
                                    δ, k, kinv,
                                    1, 0, 1, 0,
                                    0.0, 0.0)
    end
end
# Other constructors tbd

"""
    clusters(td::TDigest)
    
Extract a view as a sorted vector of the clusters contained in `td`.
"""
clusters(td::TDigest) = view(td.clusters, td._c_first:td._c_last)

"""
    buffer(td::TDigest)
    
Extract a view of the clusters buffered to be merged in `td`.
"""
buffer(td::TDigest) = view(td.buffer, td._b_first:td._b_last)

# Overloading methods for `AbstractVector`.
Base.length(td::TDigest) = td._c_last - td._c_first + 1
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

"""
    append!(td::TDigest, newdata)

Append the items of `newdata` to the t-digest `td`,
using the progressive merging algorithm.
"""
function Base.append!(td::TDigest, newdata)
    newdata = Iterators.Stateful(newdata)
    forward   = true
    newdigest = isempty(td)
    td._b_first = 1
    while length(newdata) > 0
        n = length(td)
        td.buffer[1:n] .= clusters(td) # fill the buffer with existing clusters
        insertbuffer!(td, newdata, n)  # fill with yet unused items of newdata
        bf = buffer(td)
        sort!(bf)
        # update minimum and maximum values
        if newdigest
            td._min = centroid(bf[1])
            td._max = centroid(bf[end])
            newdigest = false
        else
            (centroid(bf[1]) < td._min) && (td._min = centroid(bf[1]))
            (centroid(bf[end]) > td._max) && (td._max = centroid(bf[end]))
        end
        # merge the buffer alternating the direction
        if forward
            mergebuffer_forward!(td)

        else
            mergebuffer_backwards!(td)
        end
        forward = !forward
    end
    td._b_last = 0
    return td
end

"""
    insertbuffer!(td::TDigest, n, newdata)
    
Insert into the buffer of `td` as many items from `newdata` as possible, starting
after the index `n` of the buffer. The function stops when the buffer is full or
when the end of `newdata` is reached.
"""
function insertbuffer!(td::TDigest, newdata, n)
    d = min(length(td.buffer)-n, length(newdata))
    # td.buffer[n+1:n+d] .= Cluster.(Iterators.take(newdata, d))
    for i=1:d
        td.buffer[n+i] = (centroid=Iterators.iterate(newdata, 1)[1], count=1)
    end
    td._b_last = n+d
    return nothing
end

    
_qlimit_forward(td::TDigest, q0) = td.kinv(td.k(q0, td.δ) + 1.0, td.δ)
_qlimit_backwards(td::TDigest, q0) = td.kinv(td.k(q0, td.δ) - 1.0, td.δ)

function mergebuffer_forward!(td::TDigest)
    bf = buffer(td)
    s = sum(count, bf)
    q0 = 0.0
    qlimit = _qlimit_forward(td, q0)
    σ = bf[1]
    n = 1
    for (i,x) ∈ enumerate(view(bf, 2:length(bf)))
        q = q0 + (σ.count + x.count)/s
        if q ≤ qlimit
            σ = σ ∪ x
        else
            td.clusters[n] = σ
            if n == length(td.clusters) # interrupt and merge backwards
                # td.buffer[1:n] .= td.clusters
                # td.buffer[n+1:n+length(bf)-i] .= view(bf, i+1:length(bf))
                td.buffer[i-n+1:i] .= td.clusters
                td._b_first = i-n+1
                td._b_last  = length(bf)
                mergebuffer_backwards!(td)
                return nothing
            end
            n += 1
            q0 += σ.count/s
            qlimit = _qlimit_forward(td, q0)
            σ = x
        end
    end
    td.clusters[n] = σ
    td._c_first = 1
    td._c_last  = n
    td._b_first = 1
    td._b_last  = 0
    return nothing
end

function mergebuffer_backwards!(td::TDigest)
    bf = buffer(td)
    s = sum(count, bf)
    q0 = 1.0
    qlimit = _qlimit_backwards(td, q0)
    σ = bf[end]
    n = td.δ
    for (i,x) ∈ enumerate(view(bf, length(bf)-1:-1:1))
        q = q0 - (σ.count + x.count)/s
        if q ≥ qlimit
            σ = σ ∪ x
        else
            td.clusters[n] = σ
            if n == 1 # interrupt and merge forward
                # td.buffer[1:length(bf)-i] .= view(bf, 1:length(bf)-i)
                td.buffer[length(bf)-i+1 : length(bf)-i+td.δ] .= td.clusters
                td._b_first = 1
                td._b_last  =td.δ+length(bf)-i
                mergebuffer_forward!(td)
                return nothing
            end
            n -= 1
            q0 -= σ.count/s
            qlimit = _qlimit_backwards(td, q0)
            σ = x
        end
    end
    td.clusters[n] = σ
    td._c_first = n
    td._c_last  = td.δ
    td._b_first = 1
    td._b_last  = 0
    return nothing
end



function quantile(td::TDigest, φ::Real)
    (φ ≈ 0.0) && return td._min
    (φ ≈ 1.0) && return td._max
    nc = length(td)
    cc = clusters(td)
    ranks = cumsum(count.(cc))
    φn = φ*ranks[end]
    k = searchsortedfirst(ranks, φn)
    if k == 1      # hit first cluster
        return interpolate_centroid_first(φn, td._min, cc[1], cc[2])
    elseif k ≥ nc  # hit last cluster
        return interpolate_centroid_last(φn, cc[nc-1], cc[nc], td._max, ranks[end])
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
        # return interpolate_centroid(φn+1, w+1, Cluster(minvalue), c, c₊)
        return interpolate_centroid(φn+1, w+1, (centroid=minvalue, count=1), c, c₊)
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
