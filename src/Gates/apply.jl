### The Cliffords

function apply(gate::CliffordGate, operator, coefficient=1.0) # TODO: write tests for this
    map_array = default_clifford_map[gate.symbol]
    return applywithmap(gate, operator, coefficient, map_array)
end

function applywithmap(gate, operator, coefficient, map_array)
    operator = copy(operator)
    qinds = gate.qinds

    lookup_op = typeof(operator)(0)
    for ii in eachelement(qinds)
        lookup_op = setelement!(lookup_op, ii, getelement(operator, qinds[ii]))
    end
    sign, new_op = map_array[lookup_op]
    for ii in eachelement(qinds)
        operator = setelement!(operator, qinds[ii], getelement(new_op, ii))
    end
    coefficient = _multiplysign!(coefficient, sign)
    return operator, coefficient
end


# function _singleapply!(gate::CliffordGate, operator, coefficient)
#     local_operator = inttosymbol(getelement(operator, gate.qind))

#     relations_function = symbol_function_map[local_operator]
#     sign, new_symbol = relations_function[gate.symbol]

#     operator = setelement!(operator, gate.qinds[1], new_symbol)
#     coefficient = _multiplysign!(coefficient, sign)
#     return operator, coefficient
# end

# function _twoapply!(gate::CliffordGate, operator, coefficient)
#     qind1, qind2 = gate.qind
#     symb1 = inttosymbol(getelement(operator, qind1))
#     symb2 = inttosymbol(getelement(operator, qind2))

#     relations_function = clifford_function_map[gate.symbol]

#     sign, new_symbol1, new_symbol2 = relations_function[symb1, symb2]

#     operator = setelement!(operator, qind1, new_symbol1)
#     operator = setelement!(operator, qind2, new_symbol2)
#     coefficient = _multiplysign!(coefficient, sign)
#     return operator, coefficient
# end

function _multiplysign!(coefficient::Number, sign)
    return coefficient * sign
end

function _multiplysign!(coefficient::NumericPathProperties, sign)
    coefficient.coeff *= sign
    return coefficient
end

### The Pauli Gates  

function apply(gate::PauliGateUnion, operator, theta, coefficient=1.0)
    if commutes(gate, operator)
        return operator, coefficient
    else
        return applynoncummuting(gate, operator, theta, coefficient)
    end
end

function applynoncummuting(gate::PauliGateUnion, operator, theta, coefficient=1.0)
    coeff1 = applycos(coefficient, theta)
    sign, new_oper = getnewoperator(gate, operator)
    coeff2 = applysin(coefficient, theta; sign=sign)

    return operator, coeff1, new_oper, coeff2   # TODO: when does this ever not allocate memory?
end


function applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

function applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

function applysin(path_properties::PathProperties, theta; sign=1, kwargs...)
    # path_properties = copy(path_properties) # copy not necesasry. Was done in applycos.
    path_properties.nsins += 1
    path_properties.freq += 1

    path_properties.coeff = applysin(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

function applycos(path_properties::PathProperties, theta; sign=1, kwargs...)
    path_properties = copy(path_properties)
    path_properties.ncos += 1
    path_properties.freq += 1

    path_properties.coeff = applycos(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

function applyidentity(coeff::Number)
    return coeff
end
function applyidentity(path_properties::PathProperties)
    return path_properties
end