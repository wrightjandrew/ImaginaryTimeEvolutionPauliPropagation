abstract type Gate end

struct StaticGate <: Gate
    symbol::Symbol
    qinds::Vector{Int}
end

struct PauliGate <: Gate
    symbols::Vector{Symbol}
    qinds::Vector{Int}
end

struct FastPauliGate{T} <: Gate where {T<:Integer}
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    bitoperator::T
end

PauliGateUnion = Union{PauliGate,FastPauliGate}

tofastgates(gate::Gate, nq::Integer) = gate

function tofastgates(pauli_gate::PauliGate, nq::Integer)
    base_ops = [:I for _ in 1:nq]
    for (qind, op) in zip(pauli_gate.qinds, pauli_gate.symbols)
        base_ops[qind] = op
    end
    return FastPauliGate(pauli_gate.symbols, pauli_gate.qinds, symboltoint(base_ops))
end

function tofastgates(circ::AbstractVector)
    nq = 1
    for gate in circ
        nq = max(nq, maximum(gate.qinds))
    end
    fast_circ = Gate[]

    for gate in circ
        push!(fast_circ, tofastgates(gate, nq))
    end
    return fast_circ
end

## Cliffords  # TODO: Implement good Bitoperation Cliffords.

const H_relations = Dict(
    :I => (1, :I),
    :X => (1, :Z),
    :Y => (-1, :Y),
    :Z => (1, :X),
)

const CNOT_relations = Dict(
    (:I, :I) => (1, :I, :I),
    (:I, :X) => (1, :I, :X),
    (:I, :Y) => (1, :Z, :Y),
    (:I, :Z) => (1, :Z, :Z),
    (:X, :I) => (1, :X, :X),
    (:X, :X) => (1, :X, :I),
    (:X, :Y) => (1, :Y, :Z),
    (:X, :Z) => (-1, :Y, :Y),
    (:Y, :I) => (1, :Y, :X),
    (:Y, :X) => (1, :Y, :I),
    (:Y, :Y) => (-1, :X, :Z),
    (:Y, :Z) => (1, :X, :Y),
    (:Z, :I) => (1, :Z, :I),
    (:Z, :X) => (1, :Z, :X),
    (:Z, :Y) => (1, :I, :Y),
    (:Z, :Z) => (1, :I, :Z),
)

const ZZpihalf_relations = Dict(  # with transpose || with permutedims
    (:I, :I) => (1, :I, :I),
    (:I, :X) => (1, :Z, :Y),
    (:I, :Y) => (-1, :Z, :X),
    (:I, :Z) => (1, :I, :Z),
    (:X, :I) => (1, :Y, :Z),
    (:X, :X) => (1, :X, :X),
    (:X, :Y) => (1, :X, :Y),
    (:X, :Z) => (1, :Y, :I),
    (:Y, :I) => (-1, :X, :Z),
    (:Y, :X) => (1, :Y, :X),
    (:Y, :Y) => (1, :Y, :Y),
    (:Y, :Z) => (-1, :X, :I),
    (:Z, :I) => (1, :Z, :I),
    (:Z, :X) => (1, :I, :Y),
    (:Z, :Y) => (-1, :I, :X),
    (:Z, :Z) => (1, :Z, :Z),
)

const clifford_function_map = Dict(
    :X => H_relations,
    :CNOT => CNOT_relations,
    :ZZpihalf => ZZpihalf_relations,
)