module SingleCellProjectionsExt


export getlocalmaxima, readstereoseqbinned


using Base.Broadcast: @__dot__
using DataFrames
using SingleCellProjections: DataMatrix
using SparseArrays: sparse

using StereoSSAM


function StereoSSAM.getlocalmaxima(
    T::Type{DataMatrix},
    counts,
    localmax,
    kernel;
    genes = nothing,
)
    mat, genes, x_y, (x, y) =
        StereoSSAM._getlocalmaxima(counts, localmax, kernel; genes = genes)
    DataMatrix(
        transpose(mat),
        DataFrame(gene = genes),
        DataFrame(coord = x_y, x = x, y = y),
    )
end

function StereoSSAM.readstereoseqbinned(file, s::Integer)
    df = StereoSSAM.loadstereoseqfile(file)

    transform!(df, @. [:x, :y] => (x -> div(x - 1, s) + 1) => [:x, :y])

    cat_coord, (x, y) = StereoSSAM.categoricalcoordinates(df.x, df.y)
    select!(df, Not([:x, :y]))

    counts = sparse(df.geneID.refs, cat_coord, df.MIDCounts)

    x_y = StereoSSAM.stringcoordinates(x, y)

    DataMatrix(
        counts,
        DataFrame(gene = df.geneID.pool),
        DataFrame(bin_id = x_y, x = x, y = y),
    )
end

end # module SingleCellProjectionsExt
