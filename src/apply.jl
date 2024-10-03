
function apply(gate::StaticGate, old_operator_dict, new_operator_dict, thetas, param_ind, args...; max_weight::Real=Inf, kwargs...)
    if length(gate.qind) == 1
        _func! = _singleapply!
        checkweight = false
    else
        _func! = _twoapply!
        checkweight = true
    end
    # TODO: should we even do truncations here given that we do not increase complexity?

    for (oper, coeff) in old_operator_dict
        oper, coeff = _func!(gate, oper, coeff)

        new_operator_dict[oper] = coeff
    end
    empty!(old_operator_dict)
    return new_operator_dict, old_operator_dict, param_ind
end

function apply(gate::PauliGateUnion, operator_dict, new_operator_dict, thetas, param_ind, args...; customtruncationfunction=nothing, min_abs_coeff=0.0, kwargs...)
    theta = thetas[param_ind]


    for (oper, old_coeff) in operator_dict

        operator_dict, new_operator_dict = applystep!(gate, oper, theta, param_ind, old_coeff, operator_dict, new_operator_dict; kwargs...)  #  max_weight=here_max_weight, 

    end

    operator_dict, new_operator_dict = mergeandclear!(operator_dict, new_operator_dict)
    param_ind -= 1

    # small coeffcient truncation
    operator_dict = _removesmallcoefficients!(operator_dict, min_abs_coeff)

    if !isnothing(customtruncationfunction)
        customtruncationfunction(operator_dict, param_ind)  # changes in-place
    end

    return operator_dict, new_operator_dict, param_ind
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