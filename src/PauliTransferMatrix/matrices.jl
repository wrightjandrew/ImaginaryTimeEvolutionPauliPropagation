###
##
# This file contains implementations of certain unitary gates.
# Support gates are local unitaries, such as one-qubit, two-qubit Pauli gates.
# Note our convetion is to define the conjugate transpose of the unitary matrix.
# In the propagation, the circuit is ran in the reverse order, so defining the
# conjugate transpose of the unitary matrix is more convenient.
##
###

# module unitaries

using LinearAlgebra

const Idmat = [1 0; 0 1]
const Xmat = [0 1; 1 0]
const Ymat = [0 -1im; 1im 0]
const Zmat = [1 0; 0 -1]

const pauli_basis = [Idmat / sqrt(2), Xmat / sqrt(2), Ymat / sqrt(2), Zmat / sqrt(2)]

const _nqubit_pauli_matrices = Dict{Int,Vector{Matrix{ComplexF64}}}(1 => pauli_basis)


"""
    getpaulimatrices(nq::Int)

Compute the Pauli basis for `n` qubits.

Arguments
- `n::Int`: The number of qubits.

Returns
- `basis::Vector{Array{ComplexF64}}`: The Pauli basis for `nq` qubits.
"""
function getpaulibasis(nq::Int)
    if haskey(_nqubit_pauli_matrices, nq)
        return _nqubit_pauli_matrices[nq]
    else
        basis::Vector{Matrix{ComplexF64}} = _computepaulimatrices(nq)
        _nqubit_pauli_matrices[nq] = basis
        return basis
    end
end


function _computepaulimatrices(nq::Int)
    if nq == 1
        return pauli_basis
    else
        # Compute the Kronecker product iteratively for n qubits
        # TODO: This is currently very type-unstable, but it's fine if this is rarely called
        basis = pauli_basis
        for _ in 2:nq
            basis = [kron(p1, p2) for p2 in pauli_basis for p1 in basis]
        end
        return basis
    end
end


"""
    tomatrix(gate::PauliRotation, params)

Compute the conjugate transpose of the unitary matrix for the PauliRotation gate.

Arguments
- `gate::PauliRotation`: The PauliRotation gate.
- `params::Vector{Float64}`: A length-1 parameter for the gate.

Returns 
- `U::Array{ComplexF64}`: The conjugate transpose of the unitary matrix.
"""
function tomatrix(gate::PauliRotation, theta)

    nqubits = length(gate.qinds)
    pauli_basis_vec = getpaulibasis(nqubits)

    # The indices of the pauli matrix are sorted in ascending order
    sorted_indices = sortperm(gate.qinds)
    pauli = gate.symbols[sorted_indices]

    # These pauli matrices are normalized by sqrt(2)^nqubits
    pauli_mat = sqrt(2)^nqubits * pauli_basis_vec[symboltoint(pauli)+1]

    id = I(2^nqubits)

    U = cos(theta / 2) * id - 1.0im * sin(theta / 2) * pauli_mat

    return U
end


"""
    tomatrix(gate::TGate)

Compute the unitary matrix for a T gate.
"""
function tomatrix(::TGate)

    U = [[1 0]; [0 exp(1.0im * pi / 4)]]

    return U
end