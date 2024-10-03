function evalagainsinitial(res_dict, filterfunc)
    val = 0.0
    for (symbs, coeff) in pairs(res_dict)
        if filterfunc(symbs)
            continue
        else
            val += getnumcoeff(coeff)
        end
    end
    return val
end

function getnumcoeff(val::Real)
    return val
end

# function getnumcoeff(val::NumericPathProperties)
#     return val.coeff
# end

function zerofilter(res_dict)
    # return Dict(k => v for (k, v) in pairs(res_dict) if !annihilatesatzero(k))
    return Dict(k => v for (k, v) in pairs(res_dict) if !containsXorY(k))
end

evalagainstzero(res_dict) = evalagainsinitial(res_dict, containsXorY)

