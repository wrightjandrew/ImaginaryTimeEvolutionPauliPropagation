### The Cliffords

function apply(gate::StaticGate, operator, coefficient)   # TODO: revamp the Clifford gate approach
    if length(gate.qinds) == 1
        func = _singleapply!
    else
        func = _twoapply!
    end

    return fun(gate, copy(operator), coefficient)
end


function _singleapply!(gate::StaticGate, operator, coefficient)
    local_operator = inttosymbol(getelement(operator, gate.qind))

    relations_function = symbol_function_map[local_operator]
    sign, new_symbol = relations_function[gate.symbol]

    operator = setelement!(operator, gate.qinds[1], new_symbol)
    coefficient = _multiplysign!(coefficient, sign)
    return operator, coefficient
end

function _twoapply!(gate::StaticGate, operator, coefficient)
    qind1, qind2 = gate.qind
    symb1 = inttosymbol(getelement(operator, qind1))
    symb2 = inttosymbol(getelement(operator, qind2))

    relations_function = clifford_function_map[gate.symbol]

    sign, new_symbol1, new_symbol2 = relations_function[symb1, symb2]

    operator = setelement!(operator, qind1, new_symbol1)
    operator = setelement!(operator, qind2, new_symbol2)
    coefficient = _multiplysign!(coefficient, sign)
    return operator, coefficient
end

function _multiplysign!(coefficient::Number, sign)
    coefficient *= sign
    return coefficient
end

function _multiplysign!(coefficient::NumericPathProperties, sign)
    coefficient.coeff *= sign
    return coefficient
end

### The Pauli Gates  

# TODO: This should not be hard-coded for the merging_bfs function. Instead the merging_bfs function should call a more generic
#       version of the function that is then stored here.

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