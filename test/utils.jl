using StereoSSAM
using Test

@testset "stringcoordinates" begin
    x, y, z = Int[0, 1, 99_999], Int[0, 20, 99_999], Int[10, 100, 1_000]
    @test StereoSSAM.Utils.stringcoordinates(x) == ["0", "1", "99999"]
    @test StereoSSAM.Utils.stringcoordinates(x, y) == ["0_0", "1_20", "99999_99999"]
    @test StereoSSAM.Utils.stringcoordinates(x, y, z) ==
        ["0_0_10", "1_20_100", "99999_99999_1000"]
end

@testset "categoricalcoordinates" begin
    x, y = Int[0, 0, 1, 0, 1], Int[0, 1, 0, 0, 1]
    @test StereoSSAM.Utils.categoricalcoordinates(x, y) ==
        ([1, 2, 3, 1, 4], ([0, 0, 1, 1], [0, 1, 0, 1]))
end