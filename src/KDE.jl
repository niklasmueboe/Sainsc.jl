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
function gaussiankernel(σ::Real, r::Real)
    l = ceil(Int, 2r * σ + 1)
    Kernel.gaussian((σ, σ), (l, l))
end

function isinbounds(x, i, b)
    x += i
    x >= 1 && x <= b
end

function kde(counts::SparseMatrixCSC{T}, kernel::OffsetArray{S}) where {T<:Real,S<:Real}
    kde_dest = zeros(S, size(counts))
    kde!(counts, kernel, kde_dest)
    kde_dest
end

function kde!(
    counts::SparseMatrixCSC{T},
    kernel::OffsetArray{S},
    dest::Matrix{S},
) where {T<:Real,S<:Real}
    m, n = size(counts)

    rows = rowvals(counts)
    vals = convert.(S, nonzeros(counts))

    for col = 1:n
        for idx in nzrange(counts, col)
            row = rows[idx]
            val = vals[idx]

            row_offsets = filter(x -> isinbounds(row, x, m), axes(kernel, 1))
            col_offsets = filter(x -> isinbounds(col, x, n), axes(kernel, 2))

            for c in col_offsets, r in row_offsets
                i = row + r
                j = col + c
                @inbounds dest[i, j] += val * kernel[r, c]
            end
        end
    end

    dest
end

function kde(counts::AbstractArray{T}, kernel) where {T<:Real}
    imfilter(counts, kernel, Fill(zero(T)), Algorithm.FIR())
end


# Local maxima detection
function findlocalmaxima(counts, d::Integer, kernel)
    total_rna = kde(totalrna(counts), kernel)
    findlocalmax(total_rna, window = (2d + 1, 2d + 1))
end


# Celltype assignment
function chunk(counts, kernel::OffsetArray, sx = 500, sy = 500)
    m, n = size(first(counts))
    x1, x2 = extrema(axes(kernel, 1))
    y1, y2 = extrema(axes(kernel, 2))

    chunks = Dict{Block{2,Int},Any}()

    cols = Int[]
    rows = Int[]
    padcols = Tuple{Int,Int}[]
    padrows = Tuple{Int,Int}[]

    for j = 1:ceil(Int, n / sx)
        left = (j - 1) * sx + 1
        right = j * sx
        l = max(1, left + x1)
        r = min(n, right + x2)

        col_chunk = map(x -> x[:, l:r], counts)
        push!(cols, min(n, right) - max(1, left) + 1)
        push!(padcols, (max(0, left - l), min(0, right - r)))

        for i = 1:ceil(Int, m / sy)
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

    chunks, rows, cols, padrows, padcols
end

function calculatecosinesim(
    counts::KeyedArray,
    signatures::AbstractMatrix{T},
    kernel::AbstractMatrix{T},
) where {T<:Real}
    n_celltypes, n_genes = size(signatures)[1:2]
    n, m = size(first(counts))

    cosine = Array{T}(undef, (n, m, n_celltypes))
    kde_norm = zeros(T, (n, m))
    kde_gene = Array{T}(undef, (n, m))
    z = zero(T)

    for (g1, g2) in zip(eachindex(counts), axes(signatures, 2))
        kde_gene .= z
        kde!(counts[g1], kernel, kde_gene)
        @. kde_norm += kde_gene^2
        for (i, ct) in enumerate(axes(signatures, 1))
            if ct == 1
                @. cosine[:, :, i] = kde_gene * signatures[ct, g2]
            else
                @views @. cosine[:, :, i] += kde_gene * signatures[ct, g2]
            end
        end
    end

    @. kde_norm = sqrt(kde_norm)
    celltype_norm = map(norm, eachslice(signatures, dims = 1))

    cosine .= cosine ./ reshape(celltype_norm, 1, 1, length(celltype_norm))
    cosine, celltype = map((x -> dropdims(x, dims = 3)), findmax(cosine; dims = 3))
    cosine ./= kde_norm

    celltypemap = map((x -> x.I[3]), celltype)

    celltypemap, cosine
end

function assigncelltype(
    counts::KeyedArray,
    signatures::AbstractDataFrame,
    kernel;
    celltypes = nothing,
)
    if !isnothing(celltypes) && length(celltypes) != nrow(signatures)
        error("Length of 'celltypes' must match number of rows in 'signatures'")
    end

    genes_exist = names(signatures) .∈ [named_axiskeys(counts)[1]]
    if !all(genes_exist)
        @warn "Not all genes in 'signatures' exist in 'counts'. Missing genes will be skipped."
    end

    counts = counts(names(signatures)[genes_exist])
    signatures = Matrix{Float32}(signatures[:, genes_exist])
    convert.(Float32, kernel)

    chunked_counts, rows, cols, padrows, padcols = chunk(counts, kernel)

    T = Int

    cosine = BlockArray(undef_blocks, Matrix{eltype(kernel)}, rows, cols)
    celltypemap = BlockArray(undef_blocks, Matrix{T}, rows, cols)

    @threads for (i, (r1, r2)) in collect(enumerate(padrows))
        @threads for (j, (c1, c2)) in collect(enumerate(padcols))
            idx = Block(i, j)
            celltypemap_chunk, cosine_chunk =
                calculatecosinesim(pop!(chunked_counts, idx), signatures, kernel)

            celltypemap[idx] = celltypemap_chunk[1+r1:end+r2, 1+c1:end+c2]
            cosine[idx] = cosine_chunk[1+r1:end+r2, 1+c1:end+c2]
        end
    end

    cosine = collect(cosine)
    celltypemap = collect(celltypemap)
    celltypemap = compress(CategoricalArray{Union{Missing,T}}(celltypemap))
    @. celltypemap[isnan(cosine)] = missing

    if !isnothing(celltypes)
        celltypemap = recode(celltypemap, Dict(enumerate(celltypes))...)
    end

    celltypemap, cosine
end
