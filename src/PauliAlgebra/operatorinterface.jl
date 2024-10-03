function getelement(oper::AbstractArray{T}, index::Int) where {T}
    return oper[index]
end

function getelement(oper::Integer, index::Integer)
    return getbitelement(oper, index)
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
    return setbitelement!(oper, index, element)
end
