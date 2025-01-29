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

"""
    totransfermap(ptm)

Computes the Pauli lookup map from a Pauli Transfer Matrix (PTM).
The PTM should be the matrix representation of a gate in Pauli basis.
The returned lookup map is a vector of vectors like [(pstr1, coeff1), (pstr2, coeff2), ...]
"""
function totransfermap(ptm::Matrix)
    col_length = size(ptm)[1]
    nq = Int(log(4, col_length))
    lookupmap = Vector{Vector{Tuple{getinttype(nq),eltype(ptm)}}}(undef, col_length)
    for (colind, colvals) in enumerate(eachcol(ptm))
        lookupmap[colind] = [(rowind, ptm[rowind, colind]) for (rowind, val) in enumerate(colvals) if val != 0]
    end
    return lookupmap
end
