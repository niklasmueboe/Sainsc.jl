export readstereoseq


using AxisKeys: wrapdims
using CSV: read as readcsv
using DataFrames
using SparseArrays: sparse
using Unzip: unzip


function loadstereoseqfile(file)
    df = readcsv(
        file,
        DataFrame,
        delim = "\t",
        pool = true,
        types = Dict("x" => Int32, "y" => Int32, "MIDCounts" => Int16),
    )
    transform!(df, [:x, :y] .=> (x -> x .- (minimum(x) - one(eltype(x)))) .=> [:x, :y])

    df
end

function readstereoseq(file)
    df = loadstereoseqfile(file)

    rows = maximum(df.x)
    cols = maximum(df.y)

    genes, counts =
        (
            (key.geneID, sparse(subdf.x, subdf.y, subdf.MIDCounts, rows, cols)) for
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
