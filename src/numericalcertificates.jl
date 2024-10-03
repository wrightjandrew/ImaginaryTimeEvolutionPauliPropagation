function montecarlosampling(circ, oper::Union{Vector{Symbol},Integer}; is_reversed=false, kwargs...)
    if !is_reversed
        circ = reverse(circ)
    end
    is_valid = true
    for (ii, gate) in enumerate(circ)
        oper, is_valid = mcapply(gate, oper; kwargs...)
        if !is_valid
            return oper, false
        end
    end
    return oper, is_valid
end



function mcapply(gate::PauliGateUnion, oper, args...; max_freq::Real=Inf, max_weight::Real=Inf, min_abs_coeff=0.0, max_sins::Real=Inf, kwargs...)  # 

    if !commutes(gate, oper)
        if rand() > 0.5
            _, oper = getnewoperator(gate, oper)  # ignore in direction
            if countweight(oper; kwargs...) > max_weight  # want to stop before splitting
                return oper, false
            end
        end
    end
    return oper, true
end
