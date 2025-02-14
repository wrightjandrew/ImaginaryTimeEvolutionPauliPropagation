using LinearAlgebra
using Random
using Test

@testset "Unitaries for Miscillaneous Gates Tests" begin
    """Test the unitary matrices."""
    tol = 1e-12

    # Test using single-qubit PauliRotation gate
    @testset "PauliRotation Y" begin
        pauligate = PauliRotation(:Y, 1)
        theta = Random.randn()

        U = tomatrix(pauligate, theta)
        exptected_U = [
            [cos(theta / 2) -sin(theta / 2)];
            [sin(theta / 2) cos(theta / 2)]
        ]
        @test LinearAlgebra.norm(U - exptected_U) < tol
    end

    # Test using T gate
    @testset "TGate" begin
        tgate = TGate(1)
        U = tomatrix(tgate)
        expected_U = [
            [1 0];
            [0 exp(1.0im * pi / 4)]
        ]
        @test LinearAlgebra.norm(U - expected_U) < tol
    end

    # Test using U1Gate
    @testset "U1Gate" begin
        gate = U1Gate([0, 1])
        params = [0.1, 0.2, 0.3, 0.4, 0., 0.]
        expected_U = [
            [1 0 0 0 ]; 
            [0 0.9875353715596337 + 0.1492513737209447im 0.04694906782365546 + 0.01713774756136047im 0 ]; 
            [0 -0.045003598147776255 - 0.021739216056256186im 0.7950889010970834 + 0.6044300803163632im 0 ]; 
            [0 0 0 1]
        ]  # this is the expected unitary matrix for the given parameters (QuEst)
        U = tomatrix(gate, params)
        
        @test LinearAlgebra.norm(U - expected_U) < tol
    end
end

