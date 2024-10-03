
containsXorY(symbs::AbstractArray{Symbol}, args...) = :X in symbs || :Y in symbs
containsXorY(int::Integer, args...) = countxy(int) > 0
containsYorZ(int::Integer, args...) = countyz(int) > 0

function countweight(pauli_op::Vector{Symbol}, args...; kwargs...)
    return sum(op in (:X, :Y, :Z) for op in pauli_op)
end

function annihilatesatzero(sym_vec::AbstractVector{Symbol})
    return any(sy in (:X, :Y, :O) for sy in sym_vec)
end