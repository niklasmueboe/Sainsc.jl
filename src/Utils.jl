module Utils

export crop!, mask!, totalrna

using Base.Broadcast: @__dot__
using Base.Threads: @threads
using PooledArrays: PooledArray
using SparseArrays: dropzeros!, findnz, nonzeros, sparse
using Unzip: unzip

"""
    crop!(counts, slice)

Crop the counts to `slice`.
"""
function crop!(counts, slice)
    @threads for i in eachindex(counts)
        counts[i] = counts[i][slice...]
    end
end

"""
    mask!(counts, mask::AbstractMatrix{Bool})

Remove all counts for which `mask` is `false`.
"""
function mask!(counts, mask::AbstractMatrix{Bool})
    inv_mask = .!mask
    @threads for i in eachindex(counts)
        x, y, _ = findnz(counts[i])
        z = zero(eltype(counts[i]))
        nonzeros(counts[i])[inv_mask[[CartesianIndex(i) for i in zip(x, y)]]] .= z
        dropzeros!(counts[i])
    end
end

"""
    totalrna(counts)

Caclulate the totalrna as sum of all genes for each pixel.
"""
function totalrna(counts)
    n = length(counts)
    x, y, v = Vector{Vector}(undef, n), Vector{Vector}(undef, n), Vector{Vector}(undef, n)

    @threads for (i, c) in collect(enumerate(counts))
        x[i], y[i], v[i] = findnz(c)
    end
    return sparse((reduce(vcat, i) for i in (x, y, v))..., size(first(counts))...)
end

function categoricalcoordinates(x...)
    coordinates = PooledArray(collect(zip(x...)), Int32)
    return coordinates.refs, unzip(coordinates.pool)
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

end # module Utils