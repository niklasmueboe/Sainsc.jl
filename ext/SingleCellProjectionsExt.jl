module SingleCellProjectionsExt

using DataFrames
using SingleCellProjections: DataMatrix

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{DataMatrix}, counts, localmax, kernel; genes=nothing
)
    mat, genes, x_y, (x, y) = StereoSSAM._getlocalmaxima(
        counts, localmax, kernel; genes=genes
    )
    return DataMatrix(
        permutedims(mat), DataFrame(; gene=genes), DataFrame(; coord=x_y, x=x, y=y)
    )
end

function StereoSSAM.readstereoseqbinned(T::Type{DataMatrix}, file, s::Integer)
    counts, genes, x_y, (x, y) = StereoSSAM._readstereoseqbinned(file, s)

    return DataMatrix(counts, DataFrame(; gene=genes), DataFrame(; bin_id=x_y, x=x, y=y))
end

end # module SingleCellProjectionsExt
