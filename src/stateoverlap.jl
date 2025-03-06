### This file contains the functions to calculate the overlap between between backpropagated pauli strings and states or general operators.

"""
    overlapbyorthogonality(psum::PauliSum, orthogonalfunc::Function)

Overlap a `PauliSum` with a state or operator via function that returns true if a Pauli string is orthogonal and hence doesn't contribute.
An example `orthogonalfunc` is `containsXorY` which returns true if a Pauli string contains an X or Y Pauli.
If not orthogonal, then a Pauli string contributes with its coefficient.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(psum::PauliSum, orthogonalfunc::F) where {F<:Function}
    if length(psum) == 0
        return 0.0
    end

    val = zero(numcoefftype(psum))
    for (pstr, coeff) in psum
        if overlapbyorthogonality(pstr, orthogonalfunc)
            val += tonumber(coeff)
        end
    end
    return val
end

"""
    overlapbyorthogonality(pstr::PauliString, orthogonalfunc::Function)

Overlap a `PauliString` with a state or operator via function that returns true if the `PauliString` is orthogonal and hence has overlap 0.
 An example `orthogonalfunc` is `containsXorY` which returns true if the `PauliString` contains an X or Y Pauli.
If not orthogonal, then the overlap is the coefficient of the `PauliString`.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(pstr::PauliString, orthogonalfunc::F) where {F<:Function}
    return !orthogonalfunc(pstr) * tonumber(pstr.coeff)
end

"""
    overlapbyorthogonality(pstr::PauliString, orthogonalfunc::Function)

Overlap an integer Pauli string with a state or operator via function that returns true if the Pauli string is orthogonal and hence has overlap 0.
 An example `orthogonalfunc` is `containsXorY` which returns true if the Pauli string contains an X or Y Pauli.
If not orthogonal, then the overlap is 1.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(pstr::PauliStringType, orthogonalfunc::F) where {F<:Function}
    return !orthogonalfunc(pstr)
end

## For the typical |0> or |+> cases
"""
    overlapwithzero(psum) 

Calculates the overlap of a Pauli sum with the zero state |0><0|
"""
overlapwithzero(psum) = overlapbyorthogonality(psum, orthogonaltozero)

"""
    overlapwithzero(pstr) 

Calculates the overlap of a Pauli string with the zero state |0><0|
"""
orthogonaltozero(pstr) = containsXorY(pstr)

"""
    overlapwithplus(psum) 

Calculates the overlap of a Pauli sum with the plus state |+><+|
"""
overlapwithplus(psum) = overlapbyorthogonality(psum, orthogonaltoplus)

"""
    orthogonaltoplus(pstr) 

Calculates the overlap of a Pauli string with the plus state |+><+|
"""
orthogonaltoplus(pstr) = containsYorZ(pstr)

# eval against |±i> not implemented

"""
    overlapwithcomputational(psum::PauliSum, onebitinds)

Calculates the overlap of a Pauli sum with the computational basis state which has one-bits at all specified `indices` and zero-bits elsewhere.
For example, `overlapwithcomputational(psum, [1,2,4])` returns the overlap with `|1101000...>`
"""
function overlapwithcomputational(psum::PauliSum, onebitinds)
    val = zero(numcoefftype(psum))
    for (pstr, coeff) in psum
        val += tonumber(coeff) * _calcsignwithones(pstr, onebitinds)
    end
    return val
end

"""
    overlapwithcomputational(pstr::PauliString, onebitinds)

Calculates the overlap of a Pauli string with the computational basis state which has one-bits at all specified `onebitinds` and zero-bits elsewhere.
For example, `overlapwithcomputational(pstr, [1,2,4])` returns the overlap with `|1101000...>` and will be either zero or plus/minus `pstr.coeff`.
"""
function overlapwithcomputational(pstr::PauliString, onebitinds)
    return _calcsignwithones(pstr.term, onebitinds) * tonumber(pstr.coeff)
end

function _calcsignwithones(pstr::PauliStringType, onebitinds)

    # factor is zero unless pstr is entirely I and Z
    if containsXorY(pstr)
        return 0
    end

    # factor is +-1 per the parity of pstr's Z=3 values at the bit=1 indices
    return (-1)^count(i -> getpauli(pstr, i) == 3, onebitinds)
end

"""
    overlapwithmaxmixed(psum::PauliSum)

Calculates the overlap of a `PauliSum` with the maximally mixed state 1/2^n I.
"""
function overlapwithmaxmixed(psum::PauliSum{TT,CT}) where {TT,CT}
    if length(psum) == 0
        return 0.0
    end

    NumType = numcoefftype(psum)

    return get(psum.terms, identitypauli(TT), zero(NumType))
end

"""
    overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)

Calculates the overlap between two `PauliSum`s.
Important: We assume 'normalized' Pauli strings, i.e. such that `Tr[P^2] = 1` for any `n`-qubit Pauli string P.
If one Pauli sum represents e.g. a normalized quantum state, the result will need to be scaled by `2^n`.  
"""
function overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)
    if length(psum1) == 0 || length(psum2) == 0
        return 0.0
    end

    NumberType = numcoefftype(psum1)

    val = float(zero(NumberType))

    longer_psum = psum1.terms
    shorter_psum = psum2.terms

    # swap psums around if the other one is sparser
    if length(longer_psum) < length(shorter_psum)
        longer_psum, shorter_psum = shorter_psum, longer_psum
    end

    # looping over the shorter psum because we are only looking for collisions
    for pstr in keys(shorter_psum)
        val += tonumber(get(longer_psum, pstr, zero(NumberType))) * tonumber(get(shorter_psum, pstr, zero(NumberType)))
    end
    return val

end


"""
    filter(psum::PauliSum, filterfunc::Function)

Return a filtered `PauliSum` by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter(psum::PauliSum, filterfunc::F) where {F<:Function}
    op_dict = filter(psum.terms, filterfunc)
    return PauliSum(psum.nqubits, op_dict)
end

"""
    filter!(psum::PauliSum, filterfunc::Function)

Filter a `PauliSum` in-place by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter!(psum::PauliSum, filterfunc::F) where {F<:Function}
    filter!(psum.terms, filterfunc)
    return psum
end

# Lower-level filter function for Dict
function filter(psum::Dict, filterfunc::F) where {F<:Function}
    return Dict(k => v for (k, v) in psum if !filterfunc(k))
end

# Lower-level in-place filter function for Dict
function filter!(psum::Dict, filterfunc::F) where {F<:Function}
    for pstr in keys(psum)
        if filterfunc(pstr)
            delete!(psum, pstr)
        end
    end
    return psum
end


# returns a new filtered dictionary, but doesn't overlap with anything
"""
    zerofilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter(psum) = filter(psum, containsXorY)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter!(psum) = filter!(psum, containsXorY)

"""
    plusfilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter(psum) = filter(psum, containsYorZ)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter!(psum) = filter!(psum, containsYorZ)




