export crop!, mask!, totalrna

using AxisKeys: named_axiskeys
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using OrderedCollections: OrderedDict
using SparseArrays

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
        nonzeros(counts[i])[inv_mask[[CartesianIndex(i) for i in zip(x, y)]]] .= zero(eltype(counts[i]))
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

    batchsize = ceil(Int, length(counts) / Threads.nthreads())
    n_batches = ceil(Int, length(counts) / batchsize)
    batches = Iterators.partition(eachindex(counts), batchsize)

    vec = Vector{Vector{SparseVector}}(undef, n_batches)
    @threads for (i, batch) in collect(enumerate(batches))
        vec[i] = _kdestack(counts, coordinates, kernel, batch)
    end

    return sparse_hcat(Iterators.flatten(vec)...)
end

function categoricalcoordinates(x, y)
    coordinates = OrderedDict{Tuple{Int,Int},Int}()
    cat = 0
    for (i, xy) in enumerate(zip(x, y))
        if !haskey(coordinates, xy)
            cat += 1
            coordinates[xy] = cat
        end
    end

    # TODO multithread?
    cat_coordinate = [coordinates[xy] for xy in zip(x, y)]
    coordinate = unzip(keys(coordinates))

    return cat_coordinate, coordinate
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

function _getlocalmaxima(counts::KeyedArray, localmax, kernel; genes=nothing)
    arr = getkdeforcoordinates(counts, localmax, kernel; genes=genes)
    if isnothing(genes)
        genes = named_axiskeys(counts)[1]
    end

    coordinates = unzip(map((c -> c.I), localmax))
    coordinatestring = stringcoordinates(coordinates...)

    return arr, genes, coordinatestring, coordinates
end
