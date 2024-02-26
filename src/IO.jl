export readstereoseq, readstereoseqbinned

using AxisKeys: wrapdims
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using CSV: read as readcsv
using DataFrames
using SparseArrays: sparse

function loadstereoseqfile(file)
    countcol_name = ["MIDCounts", "MIDCount", "UMICount"]
    countcol_type = Dict(zip(countcol_name, Iterators.repeated(Int16)))

    types = merge(countcol_type, Dict("x" => Int32, "y" => Int32))

    df = readcsv(
        file, DataFrame; delim="\t", comment="#", pool=true, types=types, validate=false
    )

    for n in countcol_name
        if n in names(df)
            rename!(df, n => :count)
        end
    end

    transform!(df, [:x, :y] .=> (x -> x .- (minimum(x) - one(eltype(x)))) .=> [:x, :y])

    return df
end

"""
    readstereoseq(file)

Read StereoSeq `file` as vector of SparseMatrixCSC.
"""
function readstereoseq(file)
    df = loadstereoseqfile(file)

    rows = maximum(df.x)
    cols = maximum(df.y)

    n = length(df.geneID.pool)
    genes = Vector{eltype(df.geneID)}(undef, n)
    counts = Vector{SparseMatrixCSC{eltype(df.count)}}(undef, n)
    @threads for (i, (key, subdf)) in collect(enumerate(pairs(groupby(df, :geneID))))
        genes[i] = key.geneID
        counts[i] = sparse(subdf.x, subdf.y, subdf.count, rows, cols)
    end

    return wrapdims(counts, genes)
end

"""
    readstereoseqbinned(file, s::Integer)

Read StereoSeq `file` and aggregate locations with bin size `s`.
"""
function readstereoseqbinned(file, s::Integer)
    df = loadstereoseqfile(file)

    transform!(df, @. [:x, :y] => (x -> div(x - 1, s) + 1) => [:x, :y])

    cat_coord, (x, y) = categoricalcoordinates(df.x, df.y)
    select!(df, Not([:x, :y]))

    counts = sparse(df.geneID.refs, cat_coord, df.count)

    return (
        counts,
        DataFrame(; gene=df.geneID.pool),
        DataFrame(; bin_id=stringcoordinates(x, y), x=x, y=y),
    )
end
