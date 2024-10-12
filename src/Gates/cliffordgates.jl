struct CliffordGate <: Gate
    symbol::Symbol
    qinds::Vector{Int}
end

# TODO: verify that these are all correct
const _default_clifford_map = Dict(
    :H => [(1, 0x00), (1, 0x03), (-1, 0x02), (1, 0x01)],
    :X => [(1, 0x00), (1, 0x01), (-1, 0x02), (-1, 0x03)],
    :Y => [(1, 0x00), (-1, 0x01), (1, 0x02), (1, 0x03)],
    :Z => [(1, 0x00), (-1, 0x01), (-1, 0x02), (1, 0x03)],
    :S => [(1, 0x00), (-1, 0x02), (1, 0x01), (1, 0x03)],
    :CNOT => [(1, 0x00), (1, 0x05), (1, 0x06), (1, 0x03), (1, 0x04), (1, 0x01), (1, 0x02), (1, 0x07), (1, 0x0b), (1, 0x0e), (-1, 0x0d), (1, 0x08), (1, 0x0f), (-1, 0x0a), (1, 0x09), (1, 0x0c)],
    :ZZpihalf => [(1, 0x00), (1, 0x0e), (-1, 0x0d), (1, 0x03), (1, 0x0b), (1, 0x05), (1, 0x06), (1, 0x08), (-1, 0x07), (1, 0x09), (1, 0x0a), (-1, 0x04), (1, 0x0c), (1, 0x02), (-1, 0x01), (1, 0x0f)],
)

const default_clifford_map = deepcopy(_default_clifford_map)

function reset_clifford_map!()
    global default_clifford_map = deepcopy(_default_clifford_map)
    return
end


# const H_relations = Dict(
#     :I => (1, :I),
#     :X => (1, :Z),
#     :Y => (-1, :Y),
#     :Z => (1, :X),
# )

# const CNOT_relations = Dict(
#     (:I, :I) => (1, :I, :I),
#     (:I, :X) => (1, :I, :X),
#     (:I, :Y) => (1, :Z, :Y),
#     (:I, :Z) => (1, :Z, :Z),
#     (:X, :I) => (1, :X, :X),
#     (:X, :X) => (1, :X, :I),
#     (:X, :Y) => (1, :Y, :Z),
#     (:X, :Z) => (-1, :Y, :Y),
#     (:Y, :I) => (1, :Y, :X),
#     (:Y, :X) => (1, :Y, :I),
#     (:Y, :Y) => (-1, :X, :Z),
#     (:Y, :Z) => (1, :X, :Y),
#     (:Z, :I) => (1, :Z, :I),
#     (:Z, :X) => (1, :Z, :X),
#     (:Z, :Y) => (1, :I, :Y),
#     (:Z, :Z) => (1, :I, :Z),
# )

# const ZZpihalf_relations = Dict(  # with transpose || with permutedims
#     (:I, :I) => (1, :I, :I),
#     (:I, :X) => (1, :Z, :Y),
#     (:I, :Y) => (-1, :Z, :X),
#     (:I, :Z) => (1, :I, :Z),
#     (:X, :I) => (1, :Y, :Z),
#     (:X, :X) => (1, :X, :X),
#     (:X, :Y) => (1, :X, :Y),
#     (:X, :Z) => (1, :Y, :I),
#     (:Y, :I) => (-1, :X, :Z),
#     (:Y, :X) => (1, :Y, :X),
#     (:Y, :Y) => (1, :Y, :Y),
#     (:Y, :Z) => (-1, :X, :I),
#     (:Z, :I) => (1, :Z, :I),
#     (:Z, :X) => (1, :I, :Y),
#     (:Z, :Y) => (-1, :I, :X),
#     (:Z, :Z) => (1, :Z, :Z),
# )

# const clifford_function_map = Dict(
#     :X => H_relations,
#     :CNOT => CNOT_relations,
#     :ZZpihalf => ZZpihalf_relations,
# )