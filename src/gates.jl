abstract type Gate end

struct StaticGate <: Gate
    symbol::String
    qind::Union{Int,Tuple{Int,Int}}
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