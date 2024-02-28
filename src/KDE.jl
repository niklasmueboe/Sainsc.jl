export gaussiankernel, findlocalmaxima, kde, assigncelltype

using AxisKeys
using Base.Broadcast: @__dot__
using Base.Threads: @threads
using BlockArrays: Block, BlockArray, undef_blocks
using CategoricalArrays
using ImageFiltering:
    Kernel, imfilter, mapwindow, Fill, Algorithm, findlocalmaxima as findlocalmax
using LinearAlgebra: norm
using OffsetArrays: OffsetArray
using SparseArrays: SparseMatrixCSC

# KDE
"""
    gaussiankernel(σ::Real, r::Real)

Generate a gaussian kernel with bandwidth `σ` and radius `r*σ`
"""
function gaussiankernel(σ::Real, r::Real)
    l = ceil(Int, 2r * σ + 1)
    return Kernel.gaussian((σ, σ), (l, l))
end

"""
    kde(counts::AbstractArray{T}, kernel) where {T<:Real}

Calculate kernel density estimate.
"""
function kde(counts::AbstractArray{T}, kernel) where {T<:Real}
    return imfilter(counts, kernel, Fill(zero(T)), Algorithm.FIR())
end

function kde(counts::SparseMatrixCSC{T}, kernel::OffsetArray{S}) where {T<:Real,S<:Real}
    kde_dest = Array{S}(undef, size(counts))
    kde!(counts, kernel, kde_dest)
    return kde_dest
end

function kde!(
    counts::SparseMatrixCSC{T}, kernel::OffsetArray{S}, dest::Matrix{S}
) where {T<:Real,S<:Real}
    m, n = size(counts)

    rows = rowvals(counts)
    vals = convert.(S, nonzeros(counts))

    kernel_row_min, kernel_row_max = extrema(axes(kernel, 1))
    kernel_col_min, kernel_col_max = extrema(axes(kernel, 2))

    dest .= zero(S)
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
    return dest
end

# Local maxima detection
"""
    findlocalmaxima(counts, d::Integer, kernel)

Find local maxima of the totalRNA KDE.

Computes the [`kde`](@ref) of the [`totalrna`](@ref) to identify local maximas.
"""
function findlocalmaxima(counts, d::Integer, kernel)
    total_rna = kde(totalrna(counts), kernel)
    return findlocalmax(total_rna; window=(2d + 1, 2d + 1))
end

# Celltype assignment
function chunk(counts, kernel::OffsetArray, sx=500, sy=500)
    m, n = size(first(counts))
    x1, x2 = extrema(axes(kernel, 1))
    y1, y2 = extrema(axes(kernel, 2))

    chunks = Dict{Block{2,Int},Any}()

    cols = Int[]
    rows = Int[]
    padcols = Tuple{Int,Int}[]
    padrows = Tuple{Int,Int}[]

    for j in 1:cld(n, sx)
        left = (j - 1) * sx + 1
        right = j * sx
        l = max(1, left + x1)
        r = min(n, right + x2)

        col_chunk = map(x -> x[:, l:r], counts)
        push!(cols, min(n, right) - max(1, left) + 1)
        push!(padcols, (max(0, left - l), min(0, right - r)))

        for i in 1:cld(m, sy)
            bottom = (i - 1) * sy + 1
            top = i * sy
            b = max(1, bottom + y1)
            t = min(m, top + y2)

            chunks[Block(i, j)] = map(x -> x[b:t, :], col_chunk)
            if j == 1
                push!(rows, min(m, top) - max(1, bottom) + 1)
                push!(padrows, (max(0, bottom - b), min(0, top - t)))
            end
        end
    end

    return chunks, rows, cols, padrows, padcols
end

function calculatecosinesim(
    counts::KeyedArray, signatures::AbstractMatrix{T}, kernel::AbstractMatrix{T}
) where {T<:Real}
    n_celltypes = size(signatures)[1]
    n, m = size(first(counts))

    cosine = Array{T}(undef, (n, m, n_celltypes))
    kde_norm = zeros(T, (n, m))
    kde_gene = Array{T}(undef, (n, m))

    for (i, (g1, g2)) in enumerate(zip(eachindex(counts), axes(signatures, 2)))
        kde!(counts[g1], kernel, kde_gene)
        @. kde_norm += kde_gene^2

        weights = reshape(view(signatures, :, g2), 1, 1, n_celltypes)
            if i == 1
            @. cosine = kde_gene * weights
            else
            @. cosine += kde_gene * weights
        end
    end

    celltype_norm = map(norm, eachslice(signatures; dims=1))
    cosine ./= reshape(celltype_norm, 1, 1, n_celltypes)
    cosine, celltype = map((x -> dropdims(x; dims=3)), findmax(cosine; dims=3))

    @. cosine /= sqrt(kde_norm)
    @. cosine[iszero(kde_norm)] = 0

    celltypemap = map((x -> x.I[3]), celltype)

    return celltypemap, cosine
end

"""
    assigncelltype(
        counts::KeyedArray, signatures::AbstractDataFrame, kernel; celltypes=nothing
    )

Assign a celltype to each pixel.

The cosine similarity is calculated using `signatures` to assign the celltype with the 
highest similarity to each pixel.

# Arguments
- `signatures::AbstractDataFrame`: celltype signatures (celltypes x genes).
- `celltypes::Vector{AbstractString}=nothing`: celltype names.
"""
function assigncelltype(
    counts::KeyedArray, signatures::AbstractDataFrame, kernel; celltypes=nothing
)
    if !isnothing(celltypes) && length(celltypes) != nrow(signatures)
        error("Length of 'celltypes' must match number of rows in 'signatures'")
    end

    genes_exist = names(signatures) .∈ [named_axiskeys(counts)[1]]
    if !all(genes_exist)
        @warn "Not all genes in 'signatures' exist in 'counts'. " *
            "Missing genes will be skipped."
    end

    T = Int
    S = Float32

    counts = counts(names(signatures)[genes_exist])
    signatures = Matrix{S}(signatures[:, genes_exist])
    kernel = convert.(S, kernel)

    chunked_counts, rows, cols, padrows, padcols = chunk(counts, kernel)

    cosine = BlockArray(undef_blocks, Matrix{S}, rows, cols)
    celltypemap = BlockArray(undef_blocks, Matrix{T}, rows, cols)

    @threads for (i, (r1, r2)) in collect(enumerate(padrows))
        @threads for (j, (c1, c2)) in collect(enumerate(padcols))
            idx = Block(i, j)
            celltypemap_chunk, cosine_chunk = calculatecosinesim(
                pop!(chunked_counts, idx), signatures, kernel
            )

            celltypemap[idx] = celltypemap_chunk[(1 + r1):(end + r2), (1 + c1):(end + c2)]
            cosine[idx] = cosine_chunk[(1 + r1):(end + r2), (1 + c1):(end + c2)]
        end
    end

    cosine = collect(cosine)
    celltypemap = collect(celltypemap)
    celltypemap = compress(CategoricalArray{Union{Missing,T}}(celltypemap))
    @. celltypemap[iszero(cosine)] = missing

    if !isnothing(celltypes)
        celltypemap = recode(celltypemap, Dict(enumerate(celltypes))...)
    end

    return celltypemap, cosine
end
