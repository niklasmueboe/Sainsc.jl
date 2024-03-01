module AutomaticSingleCellToolboxExt

using AutomaticSingleCellToolbox: WsObj
using DataFrames: DataFrame, rename!
using SparseArrays: nzrange

using StereoSSAM

function StereoSSAM.getlocalmaxima(T::Type{WsObj}, counts, localmax, kernel; genes=nothing)
    X, genes, coordinates = getlocalmaxima(counts, localmax, kernel; genes=genes)

    coordinates[!, "cell_counts"] = vec(sum(X; dims=1))
    coordinates[!, "cell_features"] = map(i -> length(nzrange(X, i)), 1:size(X, 2))

    rename!(genes, Dict(:gene => "name"))

    return WsObj(Dict("raw_dat" => X), coordinates, genes, String[], Dict())
end

function StereoSSAM.readstereoseqbinned(T::Type{WsObj}, file, s::Integer)
    X, genes, coordinates = readstereoseqbinned(file, s)

    coordinates[!, "cell_counts"] = vec(sum(X; dims=1))
    coordinates[!, "cell_features"] = map(i -> length(nzrange(X, i)), 1:size(X, 2))

    rename!(genes, Dict(:gene => "name"))

    return WsObj(Dict("raw_dat" => X), coordinates, genes, String[], Dict())
end

end # module AutomaticSingleCellToolboxExt
