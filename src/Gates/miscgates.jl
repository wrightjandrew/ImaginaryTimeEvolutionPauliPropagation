### miscgates.jl
##
# A file for gates that don't fit into the other categories.
##
###


"""
    TGate(qind::Integer) <: StaticGate

Returns a T gate acting on qubit `qind`.
It acts on qubit `qind` like a `PauliRotation(:Z, qind)` with angle π/4.
"""
struct TGate <: StaticGate
    qind::Int
end