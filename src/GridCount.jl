module GridCount

export GridCounts, gridsize

using Base.Threads: @threads
using DataFrames: DataFrame, groupby, transform!
using SparseArrays: SparseMatrixCSC, sparse

"""
    AbstractGridCounts

A stack of genes each consisting of a 2D count array of the same size.
"""
abstract type AbstractGridCounts end

"""
    GridCounts{G,T<:Number} <: AbstractGridCounts

A stack of genes`::G` each consisting of a 2D count`::T` array of the same size.
"""
mutable struct GridCounts{G,T<:Number} <: AbstractGridCounts
    counts::Dict{G,SparseMatrixCSC{T}}
    shape::Tuple{Int,Int}
    @doc """
        GridCounts(counts::Dict)

    All values of the dict must be arrays of the same size.
    """
    function GridCounts{G,T}(counts::Dict{G,SparseMatrixCSC{T}}) where {G,T<:Number}
        shape = size(first(values(counts)))
        for (_, c) in counts
            if size(c) != shape
                throw(DimensionMismatch(counts, "all `counts` must have the same size"))
            end
        end
        return new{G,T}(counts, shape)
    end
end

function GridCounts(counts::Dict{G,SparseMatrixCSC{T}}) where {G,T<:Number}
    return GridCounts{G,T}(counts)
end

"""
    GridCounts(df::DataFrame)

Construct from a DataFrame with columns 'x' (`Integer`), 'y' (`Integer`), 'count', and 
'geneID' (`Pool`).
"""
function GridCounts(df::DataFrame)
    transform!(df, [:x, :y] .=> (x -> x .- (minimum(x) - one(eltype(x)))) .=> [:x, :y])

    rows = maximum(df.x)
    cols = maximum(df.y)

    n = length(df.geneID.pool)
    genes = Vector{eltype(df.geneID)}(undef, n)
    counts = Vector{SparseMatrixCSC{eltype(df.count)}}(undef, n)
    @threads for (i, (key, subdf)) in collect(enumerate(pairs(groupby(df, :geneID))))
        genes[i] = key.geneID
        counts[i] = sparse(subdf.x, subdf.y, subdf.count, rows, cols)
    end

    return GridCounts(Dict(zip(genes, counts)))
end

"""
    GridCounts(df::DataFrame, binsize::Integer)

Construct from a DataFrame with columns 'x' (`Real`), 'y' (`Real`), 'count', and 
'geneID' (`Pool`) binning the data by `binsize`.
"""
function GridCounts(df::DataFrame, binsize::Integer)
    transform!(df, [:x, :y] .=> (x -> x .- minimum(x)) .=> [:x, :y])
    transform!(df, @. [:x, :y] => (x -> div(x, binsize) + 1) => [:x, :y])

    return GridCounts(df)
end

Base.iterate(iter::GridCounts) = iterate(iter.counts)
Base.iterate(iter::GridCounts, state) = iterate(iter.counts, state)
Base.length(collection::GridCounts) = length(collection.counts)
Base.eltype(type::GridCounts) = eltype(type.counts)
Base.haskey(collection::GridCounts, key) = haskey(collection.counts, key)
Base.getindex(collection::GridCounts, key) = getindex(collection.counts, key)
Base.get(collection::GridCounts, key, default) = get(collection.counts, key, default)
Base.delete!(collection::GridCounts, key) = delete!(collection.counts, key)
Base.pop!(collection::GridCounts) = pop!(collection.counts)
Base.keys(iterator::GridCounts) = keys(iterator.counts)
Base.values(iterator::GridCounts) = values(iterator.counts)
Base.pairs(collection::GridCounts) = pairs(collection.counts)
Base.keytype(type::GridCounts) = keytype(type.counts)
Base.valtype(type::GridCounts) = valtype(type.counts)

function Base.setindex!(collection::GridCounts, value, key)
    if size(value) == collection.shape
        setindex!(collection.counts, value, key)
    else
        throw(DimensionMismatch(value, "all `counts` must have the same size"))
    end
end

function Base.show(io::IO, x::GridCounts)
    type = "GridCounts{" * string(keytype(x)) * ", " * string(eltype(valtype(x))) * "}"
    return print(io, type, " ", x.shape, " with ", length(x), " genes")
end
# Base.show(io::IO, m::MIME"text/plain", x::GridCounts) = print(io, x)

gridsize(counts::GridCounts) = counts.shape

end # module GridCount