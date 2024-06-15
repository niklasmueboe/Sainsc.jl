module Utils

export crop!, mask!, totalrna

import ..GridCount: GridCounts, gridsize

using Base.Broadcast: @__dot__
using Base.Threads: @threads
using PooledArrays: PooledArray
using SparseArrays: dropzeros!, findnz, nonzeros, sparse
using Unzip: unzip

"""
    crop!(counts::GridCounts, slice)

Crop each gene layer in counts by indexing with `slice`.
"""
function crop!(counts::GridCounts, slice)
    for (g, c) in counts
        counts.counts[g] = c[slice...]
    end
    counts.shape = size(first(values(counts)))
    return nothing
end

"""
    mask!(counts, mask::AbstractMatrix{Bool})

Remove all counts in each gene layer for which `mask` is `false`.
"""
function mask!(counts::GridCounts{<:Any,T}, mask::AbstractMatrix{Bool}) where {T<:Any}
    inv_mask = .!mask
    z = zero(T)
    @threads for c in collect(values(counts))
        x, y, _ = findnz(c)
        nonzeros(c)[inv_mask[[CartesianIndex(i) for i in zip(x, y)]]] .= z
        dropzeros!(c)
    end
end

"""
    totalrna(counts)

Caclulate the totalrna as sum of all genes for each pixel.
"""
function totalrna(counts::GridCounts)
    n = length(counts)
    x, y, v = Vector{Vector}(undef, n), Vector{Vector}(undef, n), Vector{Vector}(undef, n)

    @threads for (i, c) in collect(enumerate(values(counts)))
        x[i], y[i], v[i] = findnz(c)
    end
    return sparse((reduce(vcat, i) for i in (x, y, v))..., gridsize(counts)...)
end

function categoricalcoordinates(x...)
    coordinates = PooledArray(collect(zip(x...)), Int32)
    return coordinates.refs, unzip(coordinates.pool)
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

end # module Utils