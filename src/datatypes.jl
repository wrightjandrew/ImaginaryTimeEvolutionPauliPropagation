
abstract type PathProperties end

mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end

NumericPathProperties(coeff) = NumericPathProperties(coeff, 0, 0, 0)



import Base: *
function *(pth::PathProperties, val::Number)
    pth.coeff *= val
    return pth
end
import Base: copy
function copy(path_properties::PathProperties)
    return typeof(path_properties)(path_properties.coeff, path_properties.nsins, path_properties.ncos, path_properties.freq)
end

Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

Base.show(io::IO, pth::NumericPathProperties) = print(io, "NumericPathProperties($(pth.coeff), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")