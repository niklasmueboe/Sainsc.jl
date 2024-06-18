module IO

export readstereoseq, readstereoseqbinned, readGEMfile

import ..Utils: categoricalcoordinates, stringcoordinates

using ..GridCount: GridCounts

using Base.Broadcast: @__dot__
using CSV: read as readcsv
using DataFrames: DataFrame, Not, rename!, select!, transform!
using SparseArrays: sparse

"""
    readGEMfile(file)

Read `file` in GEM format as `DataFrames.DataFrame`.
"""
function readGEMfile(file)
    countcol_name = ["MIDCounts", "MIDCount", "UMICount"]
    countcol_type = Dict(zip(countcol_name, Iterators.repeated(Int16)))

    types = merge(countcol_type, Dict("x" => Int32, "y" => Int32))

    df = readcsv(
        file, DataFrame; delim="\t", comment="#", pool=true, types=types, validate=false
    )

    for n in countcol_name
        if n in names(df)
            rename!(df, n => :count)
            break
        end
    end

    return df
end

"""
    readstereoseq(file)

Read StereoSeq `file` as [`GridCounts`](@ref).
"""
function readstereoseq(file)
    df = readGEMfile(file)
    transform!(df, :geneID => (x -> map(y -> convert(String, y), x; pure=true)) => :geneID)
    return GridCounts(df)
end

"""
    readstereoseqbinned(file, binsize::Integer)

Read StereoSeq `file` and aggregate into bins.
"""
function readstereoseqbinned(file, binsize::Integer)
    df = readGEMfile(file)

    transform!(df, [:x, :y] .=> (x -> x .- minimum(x)) .=> [:x, :y])
    transform!(df, @. [:x, :y] => (x -> div(x, binsize) + 1) => [:x, :y])

    cat_coord, (x, y) = categoricalcoordinates(df.x, df.y)
    select!(df, Not([:x, :y]))

    counts = sparse(df.geneID.refs, cat_coord, df.count)

    return (
        counts,
        DataFrame(; gene=df.geneID.pool),
        DataFrame(; bin_id=stringcoordinates(x, y), x=x, y=y),
    )
end

end # module IO