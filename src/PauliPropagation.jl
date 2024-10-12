module PauliPropagation

using Base.Threads

include("datatypes.jl")
export PathProperties, NumericPathProperties


include("Gates/Gates.jl")
export
    Gate,
    PauliGate,
    FastPauliGate,
    tofastgates,
    CliffordGate,
    default_clifford_map,
    reset_clifford_map!,
    apply,
    applynoncummuting

include("circuits.jl")
export
    bricklayertopology,
    get2dtopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
    tfitrottercircuit,
    heisenbergtrottercircuit,
    su4ansatz,
    qcnnansatz,
    appendSU4!

include("./PauliAlgebra/PauliAlgebra.jl")
export
    inttosymbol,
    symboltoint,
    inttostring,
    getelement,
    setelement!,
    show,
    containsXorY,
    containsYorZ

include("truncations.jl")

include("Propagation/Propagation.jl")
export mergingbfs

include("stateoverlap.jl")
export evalagainstzero, evalagainsinitial, zerofilter

include("surrogate.jl")
export operatortopathdict, PauliGateNode, gettraceevalorder, expectation, resetnodes

end
