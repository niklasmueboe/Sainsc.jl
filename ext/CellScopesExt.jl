module CellScopesExt

using DataFrames
using CellScopes: RawCountObject

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{RawCountObject}, counts, localmax, kernel; genes=nothing
)
    mat, genes, x_y, (x, y) = StereoSSAM._getlocalmaxima(
        counts, localmax, kernel; genes=genes
    )

    return (
        RawCountObject(permutedims(mat), x_y, genes), DataFrame(; cell_name=x_y, x=x, y=y)
    )
end

function StereoSSAM.readstereoseqbinned(T::Type{RawCountObject}, file, s::Integer)
    counts, genes, x_y, (x, y) = StereoSSAM._readstereoseqbinned(file, s)

    return RawCountObject(counts, x_y, genes), DataFrame(; cell_name=x_y, x=x, y=y)
end

end # module CellScopesExt
