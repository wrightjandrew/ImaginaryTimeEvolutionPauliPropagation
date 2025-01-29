### miscgates.jl
##
# A file for gates that don't fit into the other categories.
##
###


## T gate
"""
    TGate(qind::Integer) <: StaticGate

Returns a T gate acting on qubit `qind`.
It acts on qubit `qind` like a `PauliRotation(:Z, qind)` with angle π/4.
"""
struct TGate <: StaticGate
    qind::Int
end

"""
    tomatrix(gate::TGate)

Compute the unitary matrix for a `TGate`.
"""
function tomatrix(::TGate)
    return _tgate_unitary
end

const _tgate_unitary = [[1 0]; [0 exp(1.0im * pi / 4)]]
