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


## TransferMapGate
# TODO: this should all be made immutable for performance
"""
    TransferMapGate(qinds::Vector{Int}, transfer_map::Vector{Vector{Tuple{TermType,CoeffType}}}) <: StaticGate

A non-parametrized `StaticGate` defined by a transfer map acting on the qubits `qinds`.
Transfer maps can be constructed manually or generated via `totransfermap()`.
"""
struct TransferMapGate{TT,CT} <: StaticGate
    qinds::Vector{Int}
    transfer_map::Vector{Vector{Tuple{TT,CT}}}

    function TransferMapGate(qinds::Vector{Int}, transfer_map::Vector{Vector{Tuple{TT,CT}}}) where {TT,CT}
        nq = length(qinds)
        @assert nq == Int(log(4, length(transfer_map))) "The length of `qinds` `n=$nq` does not match the length of the transfer map `$(length(transfer_map)) ≠ 2^$nq`."
        return new{TT,CT}(qinds, transfer_map)
    end
end