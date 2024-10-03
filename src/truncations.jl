"""
    Truncation function should return 'true' if a path should continue to propagate
    and 'false' if it should be truncated.
"""

function frequencytruncation(path_properties::T, max_l::Real) where {T<:PathProperties}
    if max_l == Inf
        return true
    end
    return path_properties.freq <= max_l
end

function frequencytruncation(coeff::Number, max_l::Real)
    return true
end


function maxsintruncation(path_data::PathProperties, max_sins::Real)
    return path_data.nsins > max_sins
end

function maxsintruncation(path_data::Real, max_sins::Real)
    return false
end


function min_coeff_condition(coeff::Float64, min_abs_coeff)
    return abs(coeff) < min_abs_coeff
end

function min_coeff_condition(node::NumericPathProperties, min_abs_coeff)
    return abs(node.coeff) < min_abs_coeff
end

function min_coeff_condition(coeff, min_abs_coeff)
    return false
end
