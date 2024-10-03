function getinttype(nq)
    nbits = 2 * nq
    if nbits <= 8
        inttype = UInt8
    elseif nbits <= 16
        inttype = UInt16
    elseif nbits <= 32
        inttype = UInt32
    elseif nbits <= 64
        inttype = UInt64
    elseif nbits <= 128
        inttype = UInt128
    elseif nbits <= 256    # These are super slow to hash
        inttype = UInt256
        # elseif nbits <= 512
        #     inttype = UInt512
        # elseif nbits <= 1024
        #     inttype = UInt1024
    else
        inttype = BigInt  # TODO: get larger Integer types that aren't BigInt
    end
    return inttype
end

@generated function alternatingmask(int::T) where {T<:Integer}
    n_bits = min(bitsize(int), 2_048)  # for 1024 qubits
    mask = int(0)
    for ii in 0:(n_bits-1)
        if ii % 2 == 0
            mask = _setbittoone(mask, ii)
        end
    end
    return mask
end


function setelement!(oper::AbstractArray{T}, index, element::T) where {T}
    oper[index] = element
    return oper
end

function setelement!(oper::AbstractArray{T}, index, element::Integer) where {T}
    return setelement!(oper, index, inttosymbol(element))
end

function setelement!(oper::Integer, index, element::Symbol)
    return setelement!(oper, index, symboltoint(element))
end


function setelement!(oper::Integer, index, element::Integer)
    bitindex = 2 * (index - 1)

    b1 = _readbit(element, 0)
    b2 = _readbit(element, 1)

    oper = _setbit(oper, bitindex, b1)
    oper = _setbit(oper, bitindex + 1, b2)
    return oper
end

function _setbit(oper::Integer, bitindex::Integer, bit::Integer)
    if bit == true  # set to one
        oper = _setbittoone(oper, bitindex)
    else
        oper = _setbittozero(oper, bitindex)
    end
    return oper
end

function _setbittoone(oper::Integer, bitindex::Integer)
    return oper | (typeof(oper)(1) << bitindex)
end

function _setbittozero(oper::Integer, bitindex::Integer)
    oper = ~oper  # flip all bits
    oper = _setbittoone(oper, bitindex)
    oper = ~oper # flip all bits back
    return oper
end

function getelement(oper::AbstractArray{T}, index::Int) where {T}
    return oper[index]
end

function getelement(oper::Integer, index::Integer)   # TODO: This function is kinda slow.
    bitindex = 2 * (index - 1)
    return ((oper >> bitindex) & UInt8(3))
end

function _readbit(oper::Integer, bitindex::Integer)
    # bitindex begins at 0
    return (oper >> bitindex) & UInt8(1)
end