export gaussiankernel, findlocalmaxima, kde, assigncelltype

using AxisKeys
using Base.Broadcast: @__dot__
using CategoricalArrays
using ImageFiltering: Kernel, imfilter, mapwindow, Fill, Algorithm
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
    m, n = size(counts)

    kde_acc = zeros(S, (m, n))

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
                @inbounds kde_acc[i, j] += val * kernel[r, c]
            end
        end
    end

    kde_acc
end

function kde(counts::AbstractArray{T}, kernel) where {T<:Real}
    imfilter(counts, kernel, Fill(zero(T)), Algorithm.FIR())
end


# Local maxima detection
function islocalmax(arr::Matrix, x::Integer, y::Integer, s::Integer, n::Integer, m::Integer)
    is = filter(i -> isinbounds(x, i, n), -s:s)
    js = filter(i -> isinbounds(y, i, m), -s:s)

    for j in js, i in is
        xi = x + i
        yj = y + j
        if @inbounds arr[x, y] < arr[xi, yj]
            return false
        end
    end
    return true
end

function findlocalmaxima(counts, d::Integer, kernel)
    total_rna = kde(totalrna(counts), kernel)

    n, m = size(total_rna)
    max_coordinates = CartesianIndex{2}[]

    for y = 1:m, x = 1:n
        if !iszero(total_rna[x, y]) && islocalmax(total_rna, x, y, d, n, m)
            push!(max_coordinates, CartesianIndex(x, y))
        end
    end

    max_coordinates
end


# Celltype assignment
function generatecelltypemap(celltype::AbstractArray{CartesianIndex{3}}, cosine)
    celltypemap = map((x -> x.I[3]), celltype)
    celltypemap =
        compress(CategoricalArray{Union{Missing,eltype(celltypemap)}}(celltypemap))
    @. celltypemap[isnan(cosine)] = missing

    celltypemap
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

    for (g1, g2) in zip(eachindex(counts), axes(signatures, 2))
        kde_gene = kde(counts[g1], kernel)
        @. kde_norm += kde_gene^2
        for (i, ct) in enumerate(axes(signatures, 1))
            if ct == 1
                cosine[:, :, i] .= kde_gene * signatures[ct, g2]
            else
                cosine[:, :, i] .+= kde_gene * signatures[ct, g2]
            end
        end
    end

    @. kde_norm = sqrt(kde_norm)
    celltype_norm = map(norm, eachslice(signatures, dims = 1))

    cosine .= cosine ./ reshape(celltype_norm, 1, 1, length(celltype_norm))
    cosine, celltype = map((x -> dropdims(x, dims = 3)), findmax(cosine; dims = 3))
    cosine = cosine ./ kde_norm

    celltypemap = generatecelltypemap(celltype, cosine)

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

    celltypemap, cosine = calculatecosinesim(
        counts(names(signatures)[genes_exist]),
        Matrix{Float32}(signatures[:, genes_exist]),
        convert.(Float32, kernel),
    )

    if !isnothing(celltypes)
        celltypemap = recode(celltypemap, Dict(enumerate(celltypes))...)
    end

    celltypemap, cosine
end
