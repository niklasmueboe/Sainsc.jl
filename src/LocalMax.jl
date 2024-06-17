module LocalMax

export getlocalmaxima, findlocalmaxima

import ..GridCount: GridCounts
import ..Utils: stringcoordinates
import ..KDE: kde!

using Base.Threads: @threads
using DataFrames: DataFrame
using ImageFiltering: findlocalmaxima as findlocalmax
using SparseArrays: SparseVector, sparsevec, sparse_hcat
using Unzip: unzip

"""
    findlocalmaxima(img, mindist::Integer; threshold::Real=0)

Find local maxima of the `img`.

The input should be a [`kde`](@ref StereoSSAM.KDE.kde) of the 
[`totalrna`](@ref StereoSSAM.GridCount.totalrna) or comparable.
"""
function findlocalmaxima(img, mindist::Integer; threshold::Real=0)
    localmax = findlocalmax(img; window=(2mindist + 1, 2mindist + 1))
    if threshold > 0
        localmax = localmax[(img .> threshold)[localmax]]
    end
    return localmax
end

function getkdeforcoordinates(counts, coordinates, kernel)
    function _kdestack(counts, coordinates, kernel, batch)
        T = eltype(kernel)
        v = Vector{SparseVector{T}}(undef, length(batch))
        kde_gene = Array{T}(undef, size(first(counts)))

        for (i, j) in enumerate(batch)
            kde!(kde_gene, counts[j], kernel)
            v[i] = @views sparsevec(kde_gene[coordinates])
        end
        return v
    end

    batchsize = cld(length(counts), Threads.nthreads())
    n_batches = cld(length(counts), batchsize)
    batches = Iterators.partition(eachindex(counts), batchsize)

    kde_coordinates = Vector{Vector{SparseVector}}(undef, n_batches)
    @threads for (i, batch) in collect(enumerate(batches))
        kde_coordinates[i] = _kdestack(counts, coordinates, kernel, batch)
    end

    return sparse_hcat(Iterators.flatten(kde_coordinates)...)
end

"""
    getlocalmaxima(counts, localmax, kernel; genes=nothing)

Load KDE with `kernel` for coordinates at `localmax`.

# Arguments
- `localmax`: vector of locations (indices).
- `genes=nothing`: vector of genes for which to calculate KDE. If `nothing` all genes are 
    used.
"""
function getlocalmaxima(counts, localmax, kernel; genes=nothing)
    if isnothing(genes)
        genes = Vector(keys(counts))
    else
        exist = genes .∈ [keys(counts)]
        if !all(exist)
            @warn "Not all `genes` exist in 'counts'. Missing `genes` will be skipped."
        end
        genes = genes[exist]
    end

    mat = getkdeforcoordinates([counts[g] for g in genes], localmax, kernel)

    x, y = unzip(Tuple.(localmax))

    return (
        permutedims(mat),
        DataFrame(; gene=genes),
        DataFrame(; id=stringcoordinates(x, y), x=x, y=y),
    )
end

end # module LocalMax