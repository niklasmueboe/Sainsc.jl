export gaussiankernel, findlocalmaxima, kde, assigncelltype

using Base.Broadcast: @__dot__
using Base.Threads: @threads
using BlockArrays: Block, BlockArray, undef_blocks
using CategoricalArrays
using DimensionalData
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
    kde!(kde_dest, counts, kernel)
    return kde_dest
end

function kde!(
    dest::Matrix{T}, counts::SparseMatrixCSC{S}, kernel::OffsetArray{T}
) where {T<:Real,S<:Real}
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

# Local maxima detection
"""
    findlocalmaxima(img, mindist::Integer; threshold::Real=0)

Find local maxima of the `img`.

The input should be a [`kde`](@ref) of the [`totalrna`](@ref) or comparable.
"""
function findlocalmaxima(img, mindist::Integer; threshold::Real=0)
    localmax = findlocalmax(img; window=(2mindist + 1, 2mindist + 1))
    if threshold > 0
        localmax = localmax[(img .> threshold)[localmax]]
    end
    return localmax
end

# Celltype assignment
function chunk_slices(i, step, n, pad)
    bound1 = (i - 1) * step + 1
    bound2 = i * step

    slice = max(1, bound1 + pad[1]):min(n, bound2 + pad[2])
    reslice = range(1 + max(0, bound1 - slice[1]); length=min(step, n - bound1 + 1))
    return slice, reslice
end

function chunk(counts, kernel::OffsetArray, sx=500, sy=500)
    m, n = size(first(counts))
    pad_x = extrema(axes(kernel, 1))
    pad_y = extrema(axes(kernel, 2))

    chunks = Dict{Block{2,Int},Any}()

    colslices = UnitRange{Int}[]
    rowslices = UnitRange{Int}[]

    for j in 1:cld(n, sx)
        slice_x, reslice_x = chunk_slices(j, sx, n, pad_x)
        col_chunk = map(x -> x[:, slice_x], counts)
        push!(colslices, reslice_x)

        for i in 1:cld(m, sy)
            slice_y, reslice_y = chunk_slices(i, sy, m, pad_y)
            chunks[Block(i, j)] = map(x -> x[slice_y, :], col_chunk)
            if j == 1
                push!(rowslices, reslice_y)
            end
        end
    end

    return chunks, rowslices, colslices
end

function calculatecosinesim(
    counts, signatures::AbstractMatrix{T}, kernel::AbstractMatrix{T}
) where {T<:Real}
    n_celltypes = size(signatures, 1)
    n, m = size(first(counts))

    kde_norm = zeros(T, (n, m))
    kde_gene = Array{T}(undef, (n, m))

    init = true
    for (g1, g2) in zip(eachindex(counts), axes(signatures, 2))
        if counts[g1] isa AbstractSparseArray && nnz(counts[g1]) == 0
            continue
        end
        kde!(kde_gene, counts[g1], kernel)
        @. kde_norm += kde_gene^2

        weights = reshape(view(signatures, :, g2), 1, 1, n_celltypes)
        if init
            init = false
            cosine::Array{T} = kde_gene .* weights
        else
            @. cosine += kde_gene * weights
        end
    end

    # fastpath if whole chunk was "empty"
    if init
        return zeros(UInt8, (n, m)), kde_norm
    end

    celltype_norm = map(norm, eachslice(signatures; dims=1))
    cosine ./= reshape(celltype_norm, 1, 1, n_celltypes)
    cosine, celltype = map((x -> dropdims(x; dims=3)), findmax(cosine; dims=3))

    @. cosine /= sqrt(kde_norm)
    @. cosine[iszero(kde_norm)] = 0

    celltypemap = map((x -> x.I[3]), celltype)
    @. celltypemap[iszero(cosine)] = 0

    return celltypemap, cosine
end

function smallestuint(n)
    uints = (UInt8, UInt16, UInt32, UInt64)
    return uints[findfirst(x -> typemax(x) >= n, uints)]
end

"""
    function assigncelltype(
        counts::DimArray{T,1}, signatures::AbstractDataFrame, kernel; celltypes=nothing
    ) where {T<:AbstractMatrix}

Assign a celltype to each pixel.

The cosine similarity is calculated using `signatures` to assign the celltype with the 
highest similarity to each pixel.

The `eltype(kernel)` will be used for calculations and `signatures` will be cast to it.

# Arguments
- `signatures::AbstractDataFrame`: celltype signatures (celltypes x genes).
- `celltypes::Vector{AbstractString}=nothing`: celltype names.
"""
function assigncelltype(
    counts::AbstractDimArray{T,1}, signatures::AbstractDataFrame, kernel; celltypes=nothing
) where {T<:AbstractMatrix}
    if !isnothing(celltypes) && length(celltypes) != nrow(signatures)
        error("Length of 'celltypes' must match number of rows in 'signatures'")
    end

    genes_exist = names(signatures) .∈ [dims(counts, 1)]
    if !all(genes_exist)
        @warn "Not all genes in 'signatures' exist in 'counts'. " *
            "Missing genes will be skipped."
    end

    S = eltype(kernel)
    U = smallestuint(nrow(signatures))

    signatures = signatures[!, genes_exist]

    @views counts = counts[At(names(signatures))]
    signatures = Matrix{S}(signatures)

    chunked_counts, rowslices, colslices = chunk(counts, kernel)
    rows, cols = length.(rowslices), length.(colslices)

    cosine = BlockArray(undef_blocks, Matrix{S}, rows, cols)
    celltypemap = BlockArray(undef_blocks, Matrix{U}, rows, cols)

    @threads for (i, r) in collect(enumerate(rowslices))
        @threads for (j, c) in collect(enumerate(colslices))
            idx = Block(i, j)
            celltypemap_chunk, cosine_chunk = calculatecosinesim(
                pop!(chunked_counts, idx), signatures, kernel
            )
            @views celltypemap[idx] = celltypemap_chunk[r, c]
            @views cosine[idx] = cosine_chunk[r, c]
        end
    end

    cosine = collect(cosine)
    celltypemap = CategoricalMatrix{Union{Missing,U},U}(collect(celltypemap))
    recode!(celltypemap, 0 => missing)

    if !isnothing(celltypes)
        celltypemap = recode(celltypemap, Dict(enumerate(celltypes))...)
    end

    return celltypemap, cosine
end
