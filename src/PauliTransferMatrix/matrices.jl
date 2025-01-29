###
##
# This files contains functions to calculate the Pauli Transfer Matrix (PTM) of a unitary matrix.
##
###

"""
    function calculateptm(U; tol=1e-15)

Calculate the Pauli Transfer Matrix (PTM) of a unitary matrix in sparse format.
Note, by default the PTM is calculated in the Heisenberg picture, 
i.e., the PTM is that of the conjugate transpose of the unitary matrix.
Arguments
- `U::Array{}`: The unitary matrix for which the PTM is calculated.
- `tol::Float64=1e-15`: The tolerance for dropping small values in the PTM.

Returns
- `ptm::Matrix`: The PTM of the conjugate transpose of unitary matrix `U`.
"""
function calculateptm(U; tol=1e-15)
    Udag = U'
    if Udag * U ≈ U * Udag ≈ I
        # The matrix is unitary, so the PTM is real
        ET = Float64
    else
        # The matrix is not unitary, so the PTM can be complex
        ET = ComplexF64
    end


    nqubits = Int(log(2, size(Udag)[1]))

    # Write ptm as a sparse matrix
    # TODO: Some PTMs can be complex
    ptm = zeros(ET, 4^nqubits, 4^nqubits)

    # The pauli basis vector is defined to be consistent with index of the pstr
    pauli_basis_vec = getpaulibasis(nqubits)

    # We are taking udag as the actual unitary, and the PTM is defined by
    # evolving P_j and take overlap with P_i
    # PTM_{ij} = Tr(udag * P_j * udag^{\dagger} * P_i)
    # in the Pauli basis, this is always real.
    for i in 1:4^nqubits
        for j in 1:4^nqubits
            # TODO: real() is a restriction to the input unitary matrix
            val = tr(Udag * pauli_basis_vec[j] * U * pauli_basis_vec[i])

            if abs(val) < tol
                continue
            end
            ptm[i, j] = ET(val)
        end
    end

    return ptm  # return the ptm as a sparse matrix
end


## Pauli basis matrices
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
