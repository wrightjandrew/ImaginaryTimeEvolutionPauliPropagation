
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


function applystep!(gate::PauliGateUnion, oper, theta, param_idx, old_coeff, operator_dict, new_operator_dict, args...; max_freq::Real=Inf, max_weight::Real=Inf, min_abs_coeff=0.0, max_sins::Real=Inf, kwargs...)


    if commutes(gate, oper)
        # old_coeff = applyidentity(old_coeff)
        return operator_dict, new_operator_dict
    end

    # TODO: remove low level details from this level
    if !frequencytruncation(old_coeff, max_freq - 1)  # want to stop before splitting
        delete!(operator_dict, oper)
        return operator_dict, new_operator_dict
    end


    coeff1 = applycos(old_coeff, theta; param_idx=param_idx)

    # TODO: new design is to delete small values after merging. Generally make this function less complicated.
    operator_dict = _updateorremove!(oper, coeff1, operator_dict; min_abs_coeff=min_abs_coeff)

    if !maxsintruncation(old_coeff, max_sins - 1)  # -1 because we want to stop the split before it happens. TODO: adapt after pulling out truncations to top level.
        sign, new_oper = getnewoperator(gate, oper)
        coeff2 = applysin(old_coeff, theta; sign=sign, param_idx=param_idx)
        new_operator_dict = _optionallyadd!(new_oper, coeff2, new_operator_dict; max_weight=max_weight, min_abs_coeff=min_abs_coeff, kwargs...)
    end  # TODO: move this extra logic onto the level of merging_bfs

    return operator_dict, new_operator_dict
end


function _optionallyadd!(oper, coeff, new_operator_dict; max_weight::Real=Inf, min_abs_coeff=0.0, kwargs...)
    # TODO: move this to merging_bfs
    if max_weight < Inf
        if countweight(oper; kwargs...) > max_weight
            return new_operator_dict
        end
    end
    new_operator_dict[oper] = coeff
    return new_operator_dict
end

function _updateorremove!(oper, coeff, operator_dict; min_abs_coeff=0.0, kwargs...)
    # TODO: move this to merging_bfs
    operator_dict[oper] = coeff
    return operator_dict
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