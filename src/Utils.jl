export crop!, mask!, totalrna


using AxisKeys: named_axiskeys
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using OrderedCollections: OrderedDict
using SparseArrays


function crop!(counts, slice)
    @threads for i in eachindex(counts)
        counts[i] = counts[i][slice...]
    end
end

function mask!(counts, mask::AbstractMatrix)
    inv_mask = .!mask
    @threads for i in eachindex(counts)
        x, y, _ = findnz(counts[i])
        nonzeros(counts[i])[inv_mask[[CartesianIndex(i) for i in zip(x, y)]]] .= zero(eltype(counts[i]))
        dropzeros!(counts[i])
    end
end

function totalrna(counts)
    x, y, v = (reduce(vcat, i) for i in unzip(findnz(c) for c in counts))
    sparse(x, y, v, size(first(counts))...) |> collect
end

function getkdeforcoordinates(counts, coordinates, kernel; genes = nothing)
    if !isnothing(genes)
        counts = counts(genes)
    end

    arr = Matrix{Float64}(undef, (length(coordinates), length(counts)))
    @threads for (j, i) in collect(enumerate(eachindex(counts)))
        arr[:, j] .= kde(counts[i], kernel)[coordinates]
    end

    arr
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

    cat_coordinate, coordinate
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

function _getlocalmaxima(counts::KeyedArray, localmax, kernel; genes = nothing)
    arr = getkdeforcoordinates(counts, localmax, kernel; genes = genes)
    if isnothing(genes)
        genes = named_axiskeys(counts)[1]
    end

    coordinates = unzip(map((c -> c.I), localmax))
    coordinatestring = stringcoordinates(coordinates...)

    arr, genes, coordinatestring, coordinates
end
