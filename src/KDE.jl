module KDE

export gaussiankernel, kde, assigncelltype

using ..GridCount: GridCounts, gridsize

using Base.Broadcast: @__dot__
using Base.Threads: @threads
using BlockArrays: Block, BlockArray, undef_blocks
using CategoricalArrays: CategoricalMatrix, recode, recode!
using DataFrames: AbstractDataFrame, nrow
using ImageFiltering: Algorithm, Fill, Kernel, imfilter, reflect
using LinearAlgebra: dot, norm
using OffsetArrays: OffsetArray
using SparseArrays: AbstractSparseArray, SparseMatrixCSC, nnz, nonzeros, nzrange, rowvals
using Unzip: unzip

# KDE
"""
    gaussiankernel(bw::Real, r::Real)

Generate a gaussian kernel with bandwidth `bw` and radius `r * bw`
"""
function gaussiankernel(bw::Real, r::Real)
    l = ceil(Int, 2r * bw + 1)
    return Kernel.gaussian((bw, bw), (l, l))
end

"""
    kde(counts, kernel)

Calculate kernel density estimate.

# Arguments
- `kernel`: usually a centered `OffsetArrays.OffsetArray`.
"""
function kde(counts, kernel)
    return imfilter(counts, reflect(kernel), Fill(zero(eltype(counts))), Algorithm.FIR())
end

function kde(counts::SparseMatrixCSC, kernel)
    kde_dest = Array{eltype(kernel)}(undef, size(counts))
    kde!(kde_dest, counts, kernel)
    return kde_dest
end

function kde!(
    dest::AbstractMatrix{T}, counts::SparseMatrixCSC, kernel::AbstractMatrix{T}
) where {T}
    m, n = size(counts)

    rows = rowvals(counts)
    vals = convert.(T, nonzeros(counts))

    kernel_row_min, kernel_row_max = extrema(axes(kernel, 1))
    kernel_col_min, kernel_col_max = extrema(axes(kernel, 2))

    dest .= zero(T)
    for col in 1:n
        c_min = max(1 - col, kernel_col_min)
        c_max = min(n - col, kernel_col_max)
        for idx in nzrange(counts, col)
            row = rows[idx]
            val = vals[idx]

            r_min = max(1 - row, kernel_row_min)
            r_max = min(m - row, kernel_row_max)

            # @inbounds
            @views @. dest[(row + r_min):(row + r_max), (col + c_min):(col + c_max)] +=
                val * kernel[r_min:r_max, c_min:c_max]
        end
    end
    return nothing
end

# Celltype assignment
function chunk_slices(i, step, n, pad)
    start = (i - 1) * step + 1
    stop = i * step

    slice = max(1, start + pad[1]):min(n, stop + pad[2])
    unpad = range(1 + max(0, start - slice[1]); length=min(step, n - start + 1))
    return slice, unpad
end

function padinfo(counts, kernel)
    padrow = extrema(axes(kernel, 1))
    padcol = extrema(axes(kernel, 2))
    return size(first(counts)), (padrow, padcol)
end

function chunkinfo(counts, kernel, chunksize)
    function chunklengths(n, s, pad)
        lengths = Int[]
        for i in 1:cld(n, s)
            _, unpad = chunk_slices(i, s, n, pad)
            push!(lengths, length(unpad))
        end
        return lengths
    end
    srow, scol = chunksize
    (m, n), (padrow, padcol) = padinfo(counts, kernel)

    return chunklengths(m, srow, padrow), chunklengths(n, scol, padcol)
end

function getchunk(counts, kernel, i, j; chunksize=(500, 500))
    (m, n), (padrow, padcol) = padinfo(counts, kernel)
    srow, scol = chunksize

    slicerow, unpadrow = chunk_slices(i, srow, m, padrow)
    slicecol, unpadcol = chunk_slices(j, scol, n, padcol)
    chunk = map(x -> x[slicerow, slicecol], counts)

    return chunk, (unpadrow, unpadcol)
end

function getchunks(counts, kernel; chunksize=(500, 500))
    (m, n), (padrow, padcol) = padinfo(counts, kernel)
    srow, scol = chunksize

    chunks = Dict{Block{2,Int},Any}()

    colslices = UnitRange{Int}[]
    rowslices = UnitRange{Int}[]

    for j in 1:cld(n, scol)
        slicecol, unpadcol = chunk_slices(j, scol, n, padcol)
        col_chunk = map(x -> x[:, slicecol], counts)
        push!(colslices, unpadcol)

        for i in 1:cld(m, srow)
            slicerow, unpadrow = chunk_slices(i, srow, m, padrow)
            chunks[Block(i, j)] = map(x -> x[slicerow, :], col_chunk)
            if j == 1
                push!(rowslices, unpadrow)
            end
        end
    end

    return chunks, rowslices, colslices
end

function findcelltypescore(cosine, upperbound)
    max = 0
    max2 = 0
    argmax = 1
    argmax2 = 1
    for (i, val) in pairs(cosine)
        if val > max2
            if val > max
                max2 = max
                max = val
                argmax2 = argmax
                argmax = i
            else
                max2 = val
                argmax2 = i
            end
        end
    end
    score = (max - max2) / upperbound[argmax, argmax2]
    return (max, score, argmax)
end

function calculatecosinesim(
    counts, signatures, similaritycorrection, kernel, unpad; log=false
)
    # signatures must be normed row-wise (celltype-wise)
    n_celltypes = size(signatures, 1)

    T = promote_type(eltype(signatures), eltype(kernel))

    kde_norm = zeros(T, length.(unpad))
    kde_gene = Array{T}(undef, size(first(counts)))

    init = true
    for (g1, g2) in zip(eachindex(counts), axes(signatures, 2))
        if counts[g1] isa AbstractSparseArray && nnz(counts[g1]) == 0
            continue
        end
        kde!(kde_gene, counts[g1], kernel)

        if log
            @. kde_gene = log1p(kde_gene)
        end

        weights = reshape(view(signatures, :, g2), 1, 1, n_celltypes)
        kde_gene_crop = @view kde_gene[unpad...]
        @. kde_norm += kde_gene_crop^2

        if init
            init = false
            cosine::Array{T} = kde_gene_crop .* weights
        else
            @. cosine += kde_gene_crop * weights
        end
    end

    # fastpath if whole chunk was "empty"
    if init
        return zeros(UInt8, length.(unpad)), kde_norm, zeros(T, length.(unpad))
    end

    cosine, score, celltype = unzip(
        map(x -> findcelltypescore(x, similaritycorrection), eachslice(cosine; dims=(1, 2)))
    )

    empty = iszero.(kde_norm)

    @. kde_norm = sqrt(kde_norm)
    cosine ./= kde_norm
    score ./= kde_norm

    # set empty pixels to zero (to avoid nan and "random" celltype)
    cosine[empty] .= 0
    score[empty] .= 0
    celltype[empty] .= 0

    return celltype, cosine, score
end

function smallestuint(n::Integer)
    if n < 0
        throw(DomainError(n, "Must be a positive integer"))
    end
    uints = (UInt8, UInt16, UInt32, UInt64)
    return uints[findfirst(x -> typemax(x) >= n, uints)]
end

"""
    assigncelltype(counts, signatures, kernel; celltypes=nothing, log=false) -> (celltypes, cosine)

Assign a celltype to each pixel.

The cosine similarity is calculated using `signatures` to assign the celltype with the 
highest similarity to each pixel.

The `eltype(kernel)` will be used for calculations and `signatures` will be cast to it.

# Arguments
- `signatures::AbstractDataFrame`: celltype signatures (celltypes x genes).
- `celltypes::Vector{AbstractString}=nothing`: optional celltype names.
- `log::Bool`: whether to log-transform the KDE. Useful if `signatures` are calculated 
    from log-transformed gene expression.
- `lowmemory::Bool`: if false will use more memory but may have shorter runtime.
"""
function assigncelltype(
    counts,
    signatures,
    kernel;
    celltypes=nothing,
    log=false,
    lowmemory=true,
    chunksize=(500, 500),
)
    if !isnothing(celltypes) && length(celltypes) != nrow(signatures)
        error("Length of 'celltypes' must match number of rows in 'signatures'")
    end

    genes = names(signatures)
    exist = genes .∈ [keys(counts)]
    if !all(exist)
        @warn "Not all genes in 'signatures' exist in 'counts'. " *
            "Missing genes will be skipped."
        genes = genes[exist]
    end

    T = eltype(kernel)
    U = smallestuint(nrow(signatures))

    signatures = Matrix{T}(signatures[!, exist])
    signatures ./= map(norm, eachslice(signatures; dims=1))

    nsignatures = size(signatures, 1)
    sigcorrection = zeros(T, nsignatures, nsignatures)
    for (i, j) in Tuple.(CartesianIndices(sigcorrection))
        if i != j
            s = signatures[i, :] .- signatures[j, :]
            @. s[s < 0] = 0
            sigcorrection[i, j] = sqrt(dot(s, s))
        end
    end

    m, n = gridsize(counts)

    counts = [counts[g] for g in genes]

    chunklengths = chunkinfo(counts, kernel, chunksize)
    cosine = BlockArray(undef_blocks, Matrix{T}, chunklengths...)
    score = BlockArray(undef_blocks, Matrix{T}, chunklengths...)
    celltypemap = BlockArray(undef_blocks, Matrix{U}, chunklengths...)

    if lowmemory
        srow, scol = chunksize
        @threads for i in 1:cld(m, srow)
            @threads for j in 1:cld(n, scol)
                idx = Block(i, j)
                chunk, unpad = getchunk(counts, kernel, i, j; chunksize=chunksize)
                celltypemap[idx], cosine[idx], score[idx] = calculatecosinesim(
                    chunk, signatures, sigcorrection, kernel, unpad; log=log
                )
            end
        end
    else
        chunked_counts, rowslices, colslices = getchunks(
            counts, kernel; chunksize=chunksize
        )

        @threads for (i, r) in collect(enumerate(rowslices))
            @threads for (j, c) in collect(enumerate(colslices))
                idx = Block(i, j)
                chunk = pop!(chunked_counts, idx)
                celltypemap[idx], cosine[idx], score[idx] = calculatecosinesim(
                    chunk, signatures, sigcorrection, kernel, (r, c); log=log
                )
            end
        end
    end

    cosine = Matrix(cosine)
    score = Matrix(score)
    celltypemap = CategoricalMatrix{Union{Missing,U},U}(collect(celltypemap))
    recode!(celltypemap, 0 => missing)

    if !isnothing(celltypes)
        celltypemap = recode(celltypemap, Dict(enumerate(celltypes))...)
    end

    return celltypemap, cosine, score
end

end # module KDE
