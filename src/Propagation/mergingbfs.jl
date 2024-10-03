
function mergingbfs(circ, op::Union{Vector{Symbol},Integer}, thetas; kwargs...)
    val = 1.0
    etype = eltype(thetas)
    d = Dict{typeof(op),etype}(deepcopy(op) => val)
    return mergingbfs(circ, d, thetas; kwargs...)
end

function mergingbfs(circ, op, path_properties::PathProperties, thetas; kwargs...)
    d = Dict{typeof(op),typeof(path_properties)}(deepcopy(op) => path_properties)
    return mergingbfs(circ, d, thetas; kwargs...)
end


function mergingbfs(circ, d, thetas; kwargs...)
    param_ind = length(thetas)

    second_d = typeof(d)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    for gate in reverse(circ)
        d, second_d, param_ind = apply(gate, d, second_d, thetas, param_ind; kwargs...)
    end
    return d
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

function _removesmallcoefficients!(operator_dict, min_abs_coeff)
    if min_abs_coeff > 0.0
        for (oper, coeff) in operator_dict
            if min_coeff_condition(coeff, min_abs_coeff)
                delete!(operator_dict, oper)
            end
        end
    end
    return operator_dict
end


# TODO: custom merging function beyond mergewith!
function mergeandclear!(operator_dict::Dict, new_operator_dict::Dict)
    mergewith!(merge, operator_dict, new_operator_dict)
    empty!(new_operator_dict)
    return operator_dict, new_operator_dict
end

function merge(val1, val2)
    return val1 + val2
end

function merge(pth1::PathProperties, pth2::PathProperties)
    pth1.coeff = merge(pth1.coeff, pth2.coeff)
    pth1.ncos = min(pth1.ncos, pth2.ncos)
    pth1.nsins = min(pth1.nsins, pth2.nsins)
    pth1.freq = min(pth1.freq, pth2.freq)
    return pth1
end

function merge(coeff1::Number, coeff2::Number)
    return coeff1 + coeff2
end