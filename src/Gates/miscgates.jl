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
The returned unitary is returned in Schrödinger picture form. 
"""
function tomatrix(::TGate)
    return _tgate_unitary
end

const _tgate_unitary = [[1 0]; [0 exp(1.0im * pi / 4)]]


# Derived struct for U1Gate
struct U1Gate <: ParametrizedGate
    nparams::Int
    qinds::Vector{Int}

    # Inner constructor to enforce nparams = 6
    function U1Gate(qinds::Vector{Int})
        new(6, qinds)
    end
end

function tomatrix(gate::U1Gate, theta)
    nqubits = length(gate.qinds)

    if length(theta) != gate.nparams
        throw(ArgumentError("The number of input parameters $(length(theta)) 
        does not match \the number of parameters $(gate.nparams)."))
    end

    θ, α, β, δ, φ1, φ2 = theta

    U = zeros(ComplexF64, 2^nqubits, 2^nqubits)

    U[1, 1] = exp(1.0im * φ1)
    U[4, 4] = exp(1.0im * φ2)
    U[2, 2] = exp(-0.5im * (α + β - 2δ)) * cos(θ/2)
    U[2, 3] = exp(0.5im * (α - β + 2δ)) * sin(θ/2)
    U[3, 2] = - exp( - 0.5im * (α - β - 2δ)) * sin(θ/2)
    U[3, 3] = exp(0.5im * (α + β + 2δ)) * cos(θ/2)

    return U
end

## TransferMapGate
# TODO: this should all be made immutable for performance
"""
    TransferMapGate(transfer_map::Vector{Vector{Tuple{TermType,CoeffType}}}, qinds::Vector{Int}) <: StaticGate

A non-parametrized `StaticGate` defined by a transfer map acting on the qubits `qinds`.
Transfer maps can be constructed manually or generated via `totransfermap()`.
"""
struct TransferMapGate{TT,CT} <: StaticGate
    transfer_map::Vector{Vector{Tuple{TT,CT}}}
    qinds::Vector{Int}

    function TransferMapGate(transfer_map::Vector{Vector{Tuple{TT,CT}}}, qinds) where {TT,CT}
        # accept anything that can be converted to a vector of integers
        qinds = vec(collect(qinds))
        nq = length(qinds)

        @assert nq == Int(log(4, length(transfer_map))) "The length of `qinds` `n=$nq` does not match the length of the transfer map `$(length(transfer_map)) ≠ 2^$nq`."

        return new{TT,CT}(transfer_map, qinds)
    end
end