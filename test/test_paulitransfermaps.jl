using LinearAlgebra
using Random
using Test

@testset "Unitaries PTM Tests" begin
    """Test the PTM for unitary matrices."""
    tol = 1e-12

    # Test using single-qubit PauliRotation gate
    @testset "PauliRotation Y" begin
        pauligate = PauliRotation(:Y, 1)

        theta = Random.rand()
        U = tomatrix(pauligate, theta)

        ptm = calculateptm(U)
        ptm_schrodinger = calculateptm(U, heisenberg=false)

        expected_ptm = [
            [1 0 0 0];
            [0 cos(theta) 0 -sin(theta)];
            [0 0 1 0];
            [0 sin(theta) 0 cos(theta)]
        ]

        @test LinearAlgebra.norm(ptm - expected_ptm) < tol
        @test LinearAlgebra.norm(ptm_schrodinger - transpose(expected_ptm)) < tol
    end

    # Test using T gate
    @testset "TGate" begin
        tgate = TGate(1)
        matrix = tomatrix(tgate)

        ptmmap = calculateptm(matrix)
        ptm_schrodinger = calculateptm(matrix, heisenberg=false)

        expected_ptm = [
            [1 0 0 0];
            [0 1 / sqrt(2) 1 / sqrt(2) 0];
            [0 -1 / sqrt(2) 1 / sqrt(2) 0];
            [0 0 0 1]
        ]

        @test LinearAlgebra.norm(ptmmap - expected_ptm) < tol
        @test LinearAlgebra.norm(ptm_schrodinger - transpose(expected_ptm)) < tol
    end
end

@testset "Test PTM simulation" begin
    pauli_rotation = PauliRotation(:X, 1)
    θ = π / 4
    U = tomatrix(pauli_rotation, θ)
    ptm = calculateptm(U)
    transfer_map = totransfermap(ptm)
    transfer_map_gate = TransferMapGate(transfer_map, 1)

    for symb in [:I, :X, :Y, :Z]
        pstr = PauliString(1, symb, 1)
        psum1 = propagate(pauli_rotation, pstr, θ)
        psum2 = propagate(transfer_map_gate, pstr)

        # it should be zero but we are seeing numerical imprecisions
        @test PauliPropagation.norm(psum1 - psum2) < 1e-15
    end

end


@testset "Test Transfer Maps and TransferMapGates" begin
    nq = 2
    gate = PauliRotation(:X, 2)
    theta = 0.3
    # get the PTmap on qubit 1 but apply on qubit 2
    pauli_rotation = PauliRotation(:X, 1)
    ptmap = totransfermap(1, pauli_rotation, theta)

    g = TransferMapGate(ptmap, 2)
    for symb in [:I, :X, :Y, :Z]
        pstr = PauliString(nq, symb, 2)
        psum1 = propagate(g, pstr)
        psum2 = propagate(gate, pstr, theta)
        @test psum1 == psum2
    end

    # test the matrix constructors
    U = tomatrix(pauli_rotation, theta)
    ptm = calculateptm(U)
    gU = TransferMapGate(U, 2)
    gptm = TransferMapGate(ptm, 2)
    # equality check
    @test all((t1[1] == t2[1]) && (t1[2] ≈ t2[2]) for (v1, v2) in zip(g.transfer_map, gU.transfer_map) for (t1, t2) in zip(v1, v2))
    @test all((t1[1] == t2[1]) && (t1[2] ≈ t2[2]) for (v1, v2) in zip(g.transfer_map, gptm.transfer_map) for (t1, t2) in zip(v1, v2))
    @test gU.qinds == g.qinds
    @test gptm.qinds == g.qinds



    nq = 2

    circuit = [CliffordGate(:CNOT, [1, 2]), CliffordGate(:X, 1), CliffordGate(:H, 2), TGate(1), TGate(2)]
    ptmap = totransfermap(nq, circuit)
    g = TransferMapGate(ptmap, (1, 2))

    pstr = PauliString(nq, [:Y, :X], [1, 2])
    psum1 = propagate(g, pstr)
    psum2 = propagate(circuit, pstr)
    @test psum1 == psum2


    nq = 5

    circuit = Gate[]
    rxlayer!(circuit, nq)
    rzlayer!(circuit, nq)
    ryylayer!(circuit, bricklayertopology(nq))

    thetas = randn(countparameters(circuit))
    static_circuit = freeze(circuit, thetas)

    ptmap = totransfermap(nq, circuit, thetas)
    static_ptmap = totransfermap(nq, static_circuit)
    @test ptmap == static_ptmap

    g = TransferMapGate(ptmap, 1:nq)

    pstr = PauliString(nq, [rand((:X, :Y, :Z)) for _ in 1:nq], 1:nq)
    psum1 = propagate(g, pstr)
    psum2 = propagate(circuit, pstr, thetas)
    @test psum1 == psum2

end
