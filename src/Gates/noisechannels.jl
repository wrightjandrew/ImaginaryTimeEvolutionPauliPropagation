### noisechannels.jl
##
# A file for noise channels. 
# In particular Pauli noise channels and amplitude damping noise.
##
###


# Depolarzing noise channel
"""
Abstract type for parametrized noise channels
"""
abstract type ParametrizedNoiseChannel <: ParametrizedGate end

"""
Abstract type for Pauli noise, i.e., noise that is diagonal in Pauli basis
"""
abstract type PauliNoise <: ParametrizedNoiseChannel end



struct DepolarizingNoise <: PauliNoise
    qind::Int

    @doc """
        DepolarizingNoise(qind::Int)

    A depolarizing noise channel acting on the qubit at index `qind`.
    Will damp X, Y, and Z Paulis equally by a factor of `1-p`.
    """
    DepolarizingNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    DepolarizingNoise(qind::Int, p::Real)

A frozen depolarizing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X, Y, and Z Paulis equally by a factor of `1-p`.
"""
function DepolarizingNoise(qind::Int, p::Real)
    _check_noise_strength(DepolarizingNoise, p)

    return FrozenGate(DepolarizingNoise(qind), p)
end

function isdamped(::DepolarizingNoise, pauli::PauliType)
    return pauli != 0
end


struct DephasingNoise <: PauliNoise
    qind::Int

    @doc """
        DephasingNoise(qind::Int)

    A dephasing noise channel acting on the qubit at index `qind`.
    Will damp X and Y Paulis equally by a factor of `1-p`.
    """
    DephasingNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    DephasingNoise(qind::Int, p::Real)

A frozen dephasing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X and Y Paulis equally by a factor of `1-p`.
"""
function DephasingNoise(qind::Int, p::Real)
    _check_noise_strength(DephasingNoise, p)

    return FrozenGate(DephasingNoise(qind), p)
end

function isdamped(::DephasingNoise, pauli::PauliType)
    return pauli == 1 || pauli == 2
end

### Individual Pauli noise channels

struct PauliXNoise <: PauliNoise
    qind::Int

    @doc """
        PauliXNoise(qind::Int)

    A Pauli-X noise channel acting on the qubit at index `qind`.
    Will damp X Paulis by a factor of `1-p`.
    """
    PauliXNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliXNoise(qind::Int, p::Real)

A frozen Pauli-X noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X Paulis by a factor of `1-p`.
"""
function PauliXNoise(qind::Int, p::Real)
    _check_noise_strength(PauliXNoise, p)

    return FrozenGate(PauliXNoise(qind), p)
end

function isdamped(::PauliXNoise, pauli::PauliType)
    return pauli == 1
end


struct PauliYNoise <: PauliNoise
    qind::Int

    @doc """
        PauliYNoise(qind::Int)

    A Pauli-Y noise channel acting on the qubit at index `qind`.
    Will damp Y Paulis by a factor of `1-p`.
    """
    PauliYNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliYNoise(qind::Int, p::Real)

A frozen Pauli-Y noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp Y Paulis by a factor of `1-p`.
"""
function PauliYNoise(qind::Int, p::Real)
    _check_noise_strength(PauliYNoise, p)

    return FrozenGate(PauliYNoise(qind), p)
end

function isdamped(::PauliYNoise, pauli::PauliType)
    return pauli == 2
end


struct PauliZNoise <: PauliNoise
    qind::Int

    @doc """
        PauliZNoise(qind::Int)

    A Pauli-Z noise channel acting on the qubit at index `qind`.
    Will damp Z Paulis by a factor of `1-p`.
    """
    PauliZNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliZNoise(qind::Int, p::Real)

A frozen Pauli-Z noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp Z Paulis by a factor of `1-p`.
"""
function PauliZNoise(qind::Int, p::Real)
    _check_noise_strength(PauliZNoise, p)

    return FrozenGate(PauliZNoise(qind), p)
end

function isdamped(::PauliZNoise, pauli::PauliType)
    return pauli == 3
end

"""
    AmplitudeDampingNoise(qind::Int)

An amplitude damping noise channel acting on the qubit at index `qind`.
Damps X and Y Paulis by a factor of sqrt(1-gamma)
and splits Z into and gamma * I and (1-gamma) * Z component (in the transposed Heisenberg picture).
"""
struct AmplitudeDampingNoise <: ParametrizedNoiseChannel
    qind::Int
end

"""
    AmplitudeDampingNoise(qind::Int, gamma::Real)

A frozen amplitude damping noise channel acting on the qubit at index `qind` with noise strength `gamma`.
Damps X and Y Paulis, and splits Z into and I and Z component (in the transposed Heisenberg picture).
"""
function AmplitudeDampingNoise(qind::Int, gamma::Real)
    _check_noise_strength(AmplitudeDampingNoise, gamma)

    return FrozenGate(AmplitudeDampingNoise(qind), gamma)
end


function _check_noise_strength(::Type{G}, p) where {G<:ParametrizedNoiseChannel}
    if !(0 <= p <= 1)
        throw(ArgumentError("$G parameter must be between 0 and 1. Got $p."))
    end
end