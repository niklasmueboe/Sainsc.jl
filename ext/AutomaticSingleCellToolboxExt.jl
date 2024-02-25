module AutomaticSingleCellToolboxExt


using AutomaticSingleCellToolbox: WsObj
using DataFrames: DataFrame
using SparseArrays: nzrange

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{WsObj},
    counts,
    localmax,
    kernel;
    genes = nothing,
)
    X, genes, x_y, (x, y) =
        StereoSSAM._getlocalmaxima(counts, localmax, kernel; genes = genes)

    cell_counts = vec(sum(X, dims = 1))
    cell_features = sum(map(i -> length(nzrange(X, i)), 1:size(X)[2]))

    WsObj(
        Dict("raw_dat" => X),
        DataFrame(
            bin_id = x_y,
            x = x,
            y = y,
            cell_counts = cell_counts,
            cell_features = cell_features,
        ),
        DataFrame(name = genes),
        String[],
        nothing,
    )
end

function StereoSSAM.readstereoseqbinned(T::Type{WsObj}, file, s::Integer)
    X, genes, x_y, (x, y) = StereoSSAM._readstereoseqbinned(file, s)

    cell_counts = vec(sum(X, dims = 1))
    cell_features = sum(map(i -> length(nzrange(X, i)), 1:size(X)[2]))

    WsObj(
        Dict("raw_dat" => X),
        DataFrame(
            bin_id = x_y,
            x = x,
            y = y,
            cell_counts = cell_counts,
            cell_features = cell_features,
        ),
        DataFrame(name = genes),
        String[],
        nothing,
    )
end

end # module AutomaticSingleCellToolboxExt
