###
##
# This file contains functions to convert a circuit or a PTM to a transfer map.
##
###



"""
    totransfermap(nq::Integer, circuit::Vector{Gate}, thetas=nothing)

Computes the Pauli transfer map from a circuit with parameters `thetas`.
`thetas` defaults to `nothing` but is required if the circuit contains parametrized gates.
The returned lookup map is a vector of vectors like [(pstr1, coeff1), (pstr2, coeff2), ...]
"""
function totransfermap(nq::Integer, circuit, thetas=nothing)

    TermType = getinttype(nq)

    # max integer to feed into the circuit
    max_integer = 4^nq - 1

    # Do one propagation per initial Pauli string on the number of qubits (can be very expensive)
    psums = [propagate(circuit, PauliString(nq, TermType(ii), 1.0), thetas; min_abs_coeff=0) for ii in 0:max_integer]

    # Convert our transfer map style, i.e., vector of vector of tuples
    return [[(TermType(paulis), coeff) for (paulis, coeff) in psum] for psum in psums]

end

"""
    totransfermap(ptm::Matrix)

Computes the Pauli transfer map from a Pauli Transfer Matrix (PTM).
The PTM should be the matrix representation of a gate in Pauli basis.
The returned lookup map is a vector of vectors like [(pstr1, coeff1), (pstr2, coeff2), ...]
"""
function totransfermap(ptm::AbstractMatrix)
    col_length = size(ptm)[1]
    nq = Int(log(4, col_length))

    lookupmap = Vector{Vector{Tuple{getinttype(nq),eltype(ptm)}}}(undef, col_length)
    for (colind, colvals) in enumerate(eachcol(ptm))
        # the Pauli integers need to to from 0 to 3, so subtract 1
        lookupmap[colind] = [(rowind - 1, ptm[rowind, colind]) for (rowind, val) in enumerate(colvals) if val != 0]
    end
    return lookupmap
end
