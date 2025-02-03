###
##
# This files contains functions to calculate the Pauli Transfer Matrix (PTM) of a unitary matrix.
##
###

"""
    function calculateptm(mat; tol=1e-15, heisenberg=true)

Calculate the Pauli Transfer Matrix (PTM) of a matrix `mat`. 
If `mat` is the matrix exponential of a (anti-)Hermitian matrix, the PTM will be real-valued. Otherwise it can be complex.
Note, by default the PTM is calculated in the -> Heisenberg picture <-, 
i.e., the PTM is that of the conjugate transpose of the  matrix.
This can be changed via the `heisenberg::Bool` keyword argument.
Arguments
- `mat::Matrix`: The unitary matrix for which the PTM is calculated.
- `tol::Float64=1e-15`: The tolerance for dropping small values in the PTM.
- `heisenberg::Bool=true`: Whether the PTM is calculated in the Heisenberg picture. 

Returns
- `ptm::Matrix`: The PTM of the conjugate transpose of matrix `mat`.
"""
function calculateptm(mat; tol=1e-15, heisenberg=true)
    mat_dag = mat'

    nqubits = Int(log(2, size(mat_dag)[1]))

    if nqubits >= 6
        @warn (
            "Entering `calculateptm()` with $nqubits qubits. " *
            "This may soon local machines out of memory or take very long."
        )
        flush(stderr)
    end

    # Initialize the PTM as complex
    # Later we will see if the PTM is real
    # TODO: Can we determine this directly from `mat`?
    ptm = zeros(ComplexF64, 4^nqubits, 4^nqubits)

    # The pauli basis vector is defined to be consistent with index of the pstr
    pauli_basis_vec = getpaulibasis(nqubits)

    # TODO: sparse PTMs?
    for i in 1:4^nqubits
        for j in 1:4^nqubits
            # The PTM is defined by evolving P_j and taking the overlap with P_i
            # i.e., how much does each Pauli transform into outer Paulis.
            val = tr(pauli_basis_vec[i] * mat * pauli_basis_vec[j] * mat_dag)

            # truncate small values
            if abs(val) < tol
                continue
            end

            ptm[i, j] = val
        end
    end

    if heisenberg
        ptm = transpose(ptm)
    end

    # we don't want to unnecessarily return a complex PTM
    if all(imag(ptm) .≈ 0)
        ptm = real(ptm)
    end

    return ptm
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
