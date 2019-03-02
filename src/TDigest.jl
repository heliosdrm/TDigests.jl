module TDigest

export Cluster, TDStruct, mergedata!, insertdata, mergebuffer

include("scales.jl")

struct Cluster # tbd: support different AbstracFloat types
    centroid::Float64
    count::Int
end
Cluster(x::T) where {T<:Real} = Cluster(x/1, 1)

centroid(c::Cluster)  = c.centroid
count(c::Cluster) = c.count

function Base.union(c1::Cluster, c2::Cluster)
    count = c1.count + c2.count
    Cluster((c1.centroid*c1.count + c2.centroid*c2.count)/count, count)
end

for comp in (:isless, :isequal)
    @eval Base.$comp(c1::Cluster, c2::Cluster) = $comp(c1.centroid, c2.centroid)
end

Base.zero(T::Type{Cluster}) = Cluster(zero(Float64),0)


struct TDStruct{F<:Function, FI<:Function}
    centroids::Vector{Cluster}
    buffer::Vector{Cluster}
    δ::Int
    k::F
    kinv::FI
    # actual centroids: .centroids[_limits[1]:_limits[2]]
    # actual buffer: .buffer[1:_limits[3]]
    _limits::Vector{Int}
    # Inner constructors sets empty centroid and buffer vectors
    function TDStruct(δ::Int, k::Function=k1, kinv::Function=k1_inv, buffersize::Int=10δ)
        new{typeof(k),typeof(kinv)}(zeros(Cluster, δ), zeros(Cluster, δ+buffersize), δ, k, kinv, [1, 0, 0])
    end
end
function TDStruct(centroids, buffer, δ, k, kinv)
    nc = length(centroids)
    if nc > δ
        throw(ArgumentError("`length(centroids)` must be greater or equal than `δ`"))
    end
    td = TDStruct(δ, k, kinv, length(centroids)+length(buffersize))
    td.centroids .= centroids
    td.buffer    .= buffer
    td._limits  .= [1, nc, nc+length(buffer)]
end
# Other constructors tbd

# Index `TDStruct` as `AbstractVector{Centroid}`
centroids(td::TDStruct) = view(td.centroids, td._limits[1]:td._limits[2])
buffer(td::TDStruct)    = view(td.buffer, 1:td._limits[3])
Base.length(td::TDStruct) = td._limits[2] - td._limits[1] + 1
for f in (:getindex, :view, :iterate)
    @eval Base.$f(td::TDStruct, args...) = $f(centroids(td), args...)
end

leftweight(td::TDStruct, i) = sum(count, view(centroids(td), 1:i-1))
rightweight(td::TDStruct, i) = sum(count, view(centroids(td), 1:i))

function ksize(td::TDStruct, i)
    n = length(td)
    qleft = leftweight(td, i)/n
    qright = qleft + td[i].count/n
    return td.k(qright) - td.k(qleft)
end
    
function mergedata!(td::TDStruct, newdata, mergeorder)
    nitems = insertbuffer!(td, newdata)
    sort!(buffer(td))
    mergebuffer!(td, mergeorder)
    return nitems
end

function insertbuffer!(td::TDStruct, newdata)
    n = length(td)
    td.buffer[1:n] .= centroids(td)
    nitems = min(length(newdata), length(td.buffer) - n)
    @inbounds for i = 1:nitems
        n += 1
        td.buffer[n] = Cluster(newdata[i])
    end
    td._limits[3] = n
    return n
end


const MergeForward = false
const MergeReverse = true

_qlimit(td::TDStruct, q0, order::Val{MergeForward}) = td.kinv(td.k(q0, td.δ) + 1.0, td.δ)
_qlimit(td::TDStruct, q0, order::Val{MergeReverse}) = td.kinv(td.k(q0, td.δ) - 1.0, td.δ)

function mergebuffer!(td::TDStruct, mergeorder::Val{MergeForward})
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
            td.centroids[n] = σ
            n += 1
            q0 += σ.count/s
            qlimit = _qlimit(td, q0, mergeorder)
            σ = x
        end
    end
    td.centroids[n] = σ
    td._limits .= [1, n, 0]
    return n #length(td)
end

function mergebuffer!(td::TDStruct, mergeorder::Val{MergeReverse})
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
            td.centroids[n] = σ
            n -= 1
            q0 -= σ.count/s
            qlimit = _qlimit(td, q0, mergeorder)
            σ = x
        end
    end
    td.centroids[n] = σ
    td._limits .= [n, td.δ, 0]
    return n #length(td)
end



end # module
