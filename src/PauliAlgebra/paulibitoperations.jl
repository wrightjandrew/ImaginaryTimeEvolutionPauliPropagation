function countweight(oper::Integer; kwargs...)
    mask = alternatingmask(oper)
    m1 = oper & mask
    m2 = oper & (mask << 1)
    res = m1 | (m2 >> 1)
    return count_ones(res)
end


function countxy(oper::Integer; kwargs...)
    mask = alternatingmask(oper)

    op = oper ⊻ (oper >> 1)
    op = op & mask

    return count_ones(op)
end

function countyz(oper::Integer; kwargs...)
    mask = alternatingmask(oper)

    op = oper & (mask << 1)

    return count_ones(op)
end

function getnewoperator(gate::PauliGateUnion, oper)
    # TODO: make this faster and potentially into bitoperations
    new_oper = copy(oper)
    total_sign = -1
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_op = paulimult(gate_sym, getelement(new_oper, qind))
        total_sign *= sign
        new_oper = setelement!(new_oper, qind, new_partial_op)
    end
    return total_sign, new_oper
end

function getnewoperator(gate::FastPauliGate, oper)
    new_oper = copy(oper)
    total_sign = -1
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_op = paulimult(gate_sym, getelement(new_oper, qind))
        total_sign *= sign
        new_oper = setelement!(new_oper, qind, new_partial_op)
    end
    return total_sign, new_oper
end


function commutes(gate_generator::AbstractArray{T}, pauli_op::AbstractArray{T})::Bool where {T}
    return sum(!commutes(o1, o2) for (o1, o2) in zip(gate_generator, pauli_op)) % 2 == 0
end

function commutes(sym1::Symbol, sym2::Symbol)::Bool
    if sym1 == :I || sym2 == :I
        return true
    else
        return sym1 == sym2
    end
end

function commutes(sym1::Symbol, pauli_ind::Integer)::Bool
    return commutes(sym1, inttosymbol(pauli_ind))
end

function commutes(gate::PauliGateUnion, oper)
    return sum(!commutes(gate_sym, getelement(oper, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

function commutes(gate::FastPauliGate, oper::Integer)
    return fastcommutes(gate.bitoperator, oper)
end


function fastcommutes(op1::Integer, op2::Integer)

    mask0 = alternatingmask(op1)
    mask1 = mask0 << 1

    # obtain the left (then right) bits of each Pauli pair, in-place
    aBits0 = mask0 & op1
    aBits1 = mask1 & op1
    bBits0 = mask0 & op2
    bBits1 = mask1 & op2

    # shift left bits to align with right bits
    aBits1 = aBits1 >> 1
    bBits1 = bBits1 >> 1

    # sets '10' at every Pauli index where individual pairs don't commute
    flags = (aBits0 & bBits1) ⊻ (aBits1 & bBits0)

    # strings commute if parity of non-commuting pairs is even
    return (count_ones(flags) % 2) == 0
end



function paulimult(sym1::Symbol, sym2::Symbol)
    ind1 = symboltoint(sym1)
    ind2 = symboltoint(sym2)
    sign, ind3 = paulimult(ind1, ind2)
    return sign, inttosymbol(ind3)
end

function paulimult(sym::Symbol, pauli_ind::Integer)
    ind = symboltoint(sym)
    sign, ind3 = paulimult(ind, pauli_ind)
    return sign, ind3
end

function paulimult(pauli1::Integer, pauli2::Integer)
    if pauli1 == 0
        return 1, pauli2
    elseif pauli2 == 0
        return 1, pauli1
    elseif pauli1 == pauli2
        return 1, 0
    else
        new_pauli = 0
        for ii in 1:3
            if ii != pauli1 && ii != pauli2
                new_pauli = ii
                break
            end
        end
        sign = levicivita(pauli1, pauli2, new_pauli)
        return sign, new_pauli
    end
end

function paulimult(op1::AbstractArray{T}, op2::AbstractArray{T}) where {T}
    total_sign = -1
    new_op = [:I for _ in eachindex(op1)]
    for (ii, (o1, o2)) in enumerate(zip(op1, op2))
        sign, new_o = paulimult(o1, o2)
        total_sign *= sign
        new_op[ii] = new_o
    end
    return total_sign, new_op
end


const levicivita_lut = cat([0 0 0; 0 0 1; 0 -1 0],
    [0 0 -1; 0 0 0; 1 0 0],
    [0 1 0; -1 0 0; 0 0 0];
    dims=3)

function levicivita(n1::Integer, n2::Integer, n3::Integer)
    return levicivita_lut[n1, n2, n3]
end
