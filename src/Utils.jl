export crop!, mask!, totalRna


using AxisKeys: named_axiskeys
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using OrderedCollections: OrderedDict
using SparseArrays


function crop!(sp_c, slice)
    @threads for i in eachindex(sp_c)
        sp_c[i] = sp_c[i][slice...]
    end
end

function mask!(sp_c, mask::AbstractMatrix)
    inv_mask = .!mask
    @threads for i in eachindex(sp_c)
        x, y, _ = findnz(sp_c[i])
        nonzeros(sp_c[i])[inv_mask[[CartesianIndex(i) for i in zip(x, y)]]] .= zero(eltype(sp_c[i]))
        dropzeros!(sp_c[i])
    end
end

totalRna = sum


function getkdeforcoordinates(sp_c, coordinates, kernel; genes = nothing)
    if !isnothing(genes)
        sp_c = sp_c(genes)
    end

    arr = Matrix{Float64}(undef, (length(coordinates), length(sp_c)))
    @threads for (j, i) in collect(enumerate(eachindex(sp_c)))
        arr[:, j] = kde(sp_c[i], kernel)[coordinates]
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

function _getlocalmaxima(sp_c, localmax, kernel; genes = nothing)
    arr = getkdeforcoordinates(sp_c, localmax, kernel; genes = genes)
    if isnothing(genes)
        genes = named_axiskeys(sp_c)[1]
    end

    coordinates = unzip(map((c -> c.I), localmax))
    coordinatestring = stringcoordinates(coordinates...)

    arr, genes, coordinatestring, coordinates
end
