using Sainsc
using Test

using OffsetArrays: centered
using SparseArrays: sparse

input = transpose(
    [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    ],
)

result_kde = transpose(
    [
        [1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0.5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0.25, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0];;
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 0]
    ],
)

kernel = centered([
    [0.5, 0.5, 0.25];;
    [0, 1, 0];;
    [0, 0, 0]
])

@testset "smallestuint" begin
    @test Sainsc.KDE.smallestuint(0) == UInt8
    @test Sainsc.KDE.smallestuint(7) == UInt8
    @test Sainsc.KDE.smallestuint(255) == UInt8
    @test Sainsc.KDE.smallestuint(256) == UInt16
    @test Sainsc.KDE.smallestuint(100_000) == UInt32
    @test_throws DomainError Sainsc.KDE.smallestuint(-1)
end

@testset "gaussiankernel" begin
    @test size(gaussiankernel(2, 2)) == (9, 9)
    @test size(gaussiankernel(8, 2)) == (33, 33)
end

@testset "kde" begin
    @test kde(input, kernel) == result_kde
    @test kde(sparse(input), kernel) == result_kde
end

@testset "getkdeforcoordinates" begin
    kdeatcoord = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 2.0, 1.0, 2.0], 3, 2)

    counts = [sparse(input), sparse(input)]
    coords = map(CartesianIndex, [(1, 1), (2, 3), (2, 4)])

    Sainsc.LocalMax.getkdeforcoordinates(counts, coords, kernel) == kdeatcoord
end
