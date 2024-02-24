export readstereoseq


using AxisKeys: wrapdims
using CSV: read as readcsv
using DataFrames
using SparseArrays: sparse
using Unzip: unzip


function loadstereoseqfile(file)
    countcol_name = ["MIDCounts", "MIDCount", "UMICount"]
    countcol_type = Dict(zip(countcol_name, Iterators.repeated(Int16)))

    types = merge(countcol_type, Dict("x" => Int32, "y" => Int32))

    df = readcsv(
        file,
        DataFrame,
        delim = "\t",
        comment = "#",
        pool = true,
        types = types,
        validate = false,
    )

    for n in countcol_name
        if n in names(df)
            rename!(df, n => :count)
        end
    end

    transform!(df, [:x, :y] .=> (x -> x .- (minimum(x) - one(eltype(x)))) .=> [:x, :y])

    df
end

function readstereoseq(file)
    df = loadstereoseqfile(file)

    rows = maximum(df.x)
    cols = maximum(df.y)

    genes, counts =
        (
            (key.geneID, sparse(subdf.x, subdf.y, subdf.count, rows, cols)) for
            (key, subdf) in pairs(groupby(df, :geneID))
        ) |>
        unzip |>
        collect

    #     n = unique(df.geneID) |> length
    #     genes = Vector{AbstractString}(undef, n)
    #     counts = Vector{SparseMatrixCSC}(undef, n)
    #     @threads for (key, subdf) in pairs(groupby(df, :geneID))
    #         push!(genes, key.geneID)
    #         push!(counts, sparse(subdf.x, subdf.y, subdf.MIDCounts, rows, cols))
    #     end

    wrapdims(counts, genes)
end

function _readstereoseqbinned(file, s::Integer)
    df = loadstereoseqfile(file)

    transform!(df, @. [:x, :y] => (x -> div(x - 1, s) + 1) => [:x, :y])

    cat_coord, (x, y) = categoricalcoordinates(df.x, df.y)
    select!(df, Not([:x, :y]))

    counts = sparse(df.geneID.refs, cat_coord, df.count)

    x_y = stringcoordinates(x, y)

    counts, df.geneID.pool, x_y, (x, y)
end
