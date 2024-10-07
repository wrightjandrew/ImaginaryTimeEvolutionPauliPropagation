using PauliPropagation
using Test

@testset "PauliPropagation.jl" begin
    # Write your tests here.
    include("test_apply.jl")
    th = randn()
    applypauligate([:Z], [1], th) .== (0x00c0, 1.0)
    th = randn()
    applypauligate([:X], [1], th) .≈ (0x0003, cos(th), 0x0002, sin(th))
    th = randn()
    applypauligate([:Y], [1], th) .≈ (0x0003, cos(th), 0x0001, -sin(th))


    include("test_mergingbfs.jl")
    numericalPP(8, 4, Inf, 0.0) ≈ 0.217200
    hybridPP(8, 4, Inf, 0.0, Inf) ≈ 0.217200
    surrogatePP(8, 4, Inf, Inf) ≈ 0.217200
end


