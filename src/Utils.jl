export crop!, getlocalmaxima, mask!, totalrna

using AxisKeys: named_axiskeys
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using DataFrames: DataFrame
using PooledArrays
using SparseArrays
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
    x, y, v = (reduce(vcat, i) for i in unzip(findnz(c) for c in counts))
    return collect(sparse(x, y, v, size(first(counts))...))
end

function getkdeforcoordinates(counts, coordinates, kernel; genes=nothing)
    function _kdestack(counts, coordinates, kernel, batch)
        T = eltype(kernel)
        v = Vector{SparseVector{T}}(undef, length(batch))
        kde_gene = Array{T}(undef, size(first(counts)))

        for (i, gene) in enumerate(batch)
            kde!(counts[gene], kernel, kde_gene)
            v[i] = sparsevec(kde_gene[coordinates])
        end
        return v
    end

    if !isnothing(genes)
        counts = counts(genes)
    end

    batchsize = cld(length(counts), Threads.nthreads())
    n_batches = cld(length(counts), batchsize)
    batches = Iterators.partition(eachindex(counts), batchsize)

    vec = Vector{Vector{SparseVector}}(undef, n_batches)
    @threads for (i, batch) in collect(enumerate(batches))
        vec[i] = _kdestack(counts, coordinates, kernel, batch)
    end

    return sparse_hcat(Iterators.flatten(vec)...)
end

function categoricalcoordinates(x, y)
    coordinates = PooledArray(collect(zip(x, y)), Int32)
    return coordinates.refs, coordinates.pool
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

"""
    getlocalmaxima(counts, localmax, kernel; genes=nothing)

Load KDE with `kernel` for coordinates at `localmax`.

# Arguments
- `genes::Vector{AbstractString}=nothing`: vector of genes for which to calculate KDE.
"""
function getlocalmaxima(counts, localmax, kernel; genes=nothing)
    mat = getkdeforcoordinates(counts, localmax, kernel; genes=genes)
    if isnothing(genes)
        genes = named_axiskeys(counts)[1]
    end

    x, y = unzip(map((c -> c.I), localmax))

    return (
        permutedims(mat),
        DataFrame(; gene=genes),
        DataFrame(; id=stringcoordinates(x, y), x=x, y=y),
    )
end
