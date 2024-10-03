
abstract type PathProperties end

mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end


function PathProperties(coeff::Number)
    return NumericPathProperties(coeff, 0, 0, 0)
end


import Base: *
function *(pth::PathProperties, val::Number)
    pth.coeff *= val
    return pth
end

Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

# TODO: This is for the Surrogate
# mutable struct NodePathProperties <: PathProperties
#     coeff::CircuitNode
#     nsins::Int
#     ncos::Int
#     freq::Int
# end

# function PathProperties(node::CircuitNode)
#     return NodePathProperties(node, 0, 0, 0)
# end

# function Base.copy(pth::NodePathProperties)
#     return NodePathProperties(pth.coeff, pth.nsins, pth.ncos, pth.freq)
# end