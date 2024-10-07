using PauliPropagation
using Test

@testset "PauliPropagation.jl" begin
    # Write your tests here.
    include("test_mergingbfs.jl")
    numericalPP(8, 4, Inf, 0.0) ≈ 0.217200
    hybridPP(8, 4, Inf, 0.0, Inf) ≈ 0.217200
    surrogatePP(8, 4, Inf, Inf) ≈ 0.217200
end


