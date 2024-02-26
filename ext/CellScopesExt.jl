module CellScopesExt

using DataFrames: rename!
using CellScopes: RawCountObject

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{RawCountObject}, counts, localmax, kernel; genes=nothing
)
    mat, genes, coordinates = getlocalmaxima(counts, localmax, kernel; genes=genes)

    rename!(coordinates, Dict(:id => "cell_name"))

    return RawCountObject(mat, coordinates.cell_name, genes.gene), coordinates
end

function StereoSSAM.readstereoseqbinned(T::Type{RawCountObject}, file, s::Integer)
    counts, genes, coordinates = readstereoseqbinned(file, s)

    rename!(coordinates, Dict(:id => "cell_name"))

    return RawCountObject(counts, coordinates.cell_name, genes.gene), coordinates
end

end # module CellScopesExt
