export gaussiankernel, findlocalmaxima, kde, assigncelltype


using ImageFiltering: Kernel, imfilter, mapwindow, Fill
using SparseArrays: SparseMatrixCSC


function kde(arr::SparseMatrixCSC, kernel)
    kde(collect(arr), kernel)
end

function kde(arr, kernel)
    imfilter(collect(arr), kernel, Fill(zero(eltype(arr))))
end

function findlocalmaxima(arr, d::Integer)
    window = Tuple(2d + 1 for _ = 1:length(size(arr)))
    arr_max = mapwindow(maximum, arr, window, border = Fill(zero(eltype(arr)))) # TODO slow
    findall(.!iszero.(arr) .& (arr .== arr_max))
end

function assigncelltype(sp_c, signatures; kwargs...)

end

function gaussiankernel(σ::Real, r::Real)
    l = ceil(Int, 2r * σ + 1)
    Kernel.gaussian((σ, σ), (l, l))
end
