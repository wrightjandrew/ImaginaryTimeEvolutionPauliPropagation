
function mergingbfs(circ, op::Union{Vector{Symbol},Integer}, thetas; kwargs...)
    # if typeof(thetas) <: AbstractArray{Num}
    #     val = PathProperties(1.0)
    #     etype = PathProperties
    # else
    val = 1.0
    etype = eltype(thetas)
    # end
    d = Dict{typeof(op),etype}(deepcopy(op) => val)
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