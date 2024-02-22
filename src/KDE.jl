export gaussiankernel, findlocalmaxima, kde, assigncelltype


using ImageFiltering: Kernel, imfilter, mapwindow, Fill
using SparseArrays: SparseMatrixCSC


function kde(arr::SparseMatrixCSC, kernel)
    kde(collect(arr), kernel)
end

function kde(arr, kernel)
    imfilter(collect(arr), kernel, Fill(zero(eltype(arr))))
end

function findlocalmaxima(sp_c, d::Integer, kernel)
    function islocalmax(
        arr::Matrix,
        x::Integer,
        y::Integer,
        s::Integer,
        n::Integer,
        m::Integer,
    )
        for i = -s:s
            xi = x + i
            if xi < 1 || xi > n
                continue
            end
            for j = -s:s
                yj = y + j
                if yj < 1 || yj > m
                    continue
                elseif @inbounds arr[x, y] < arr[xi, yj]
                    return false
                end
            end
        end
        return true
    end

    total_rna = kde(totalrna(sp_c), kernel) |> collect

    n, m = size(total_rna)
    max_coordinates = CartesianIndex{2}[]

    for x = 1:n, y = 1:m
        if !iszero(total_rna[x, y]) && islocalmax(total_rna, x, y, d, n, m)
            push!(max_coordinates, CartesianIndex(x, y))
        end
    end

    max_coordinates
end

function assigncelltype(sp_c, signatures; kwargs...)

end

function gaussiankernel(σ::Real, r::Real)
    l = ceil(Int, 2r * σ + 1)
    Kernel.gaussian((σ, σ), (l, l))
end
