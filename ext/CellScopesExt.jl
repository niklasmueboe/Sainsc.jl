module CellScopesExt

using DataFrames: rename!
using CellScopes: RawCountObject

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    ::Type{RawCountObject}, counts, localmax, kernel; genes=nothing
)
    mat, genes, coordinates = getlocalmaxima(counts, localmax, kernel; genes=genes)

    rename!(coordinates, Dict(:id => "cell_name"))

    return RawCountObject(mat, coordinates.cell_name, genes.gene), coordinates
end

function StereoSSAM.readstereoseqbinned(::Type{RawCountObject}, file, binsize::Integer)
    counts, genes, coordinates = readstereoseqbinned(file, binsize)

    rename!(coordinates, Dict(:id => "cell_name"))

    return RawCountObject(counts, coordinates.cell_name, genes.gene), coordinates
end

end # module CellScopesExt
