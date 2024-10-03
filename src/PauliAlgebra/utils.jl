
const pauli_ops::Vector{Symbol} = [:I, :X, :Y, :Z]

symboltoint(sym::Symbol) = findfirst(s -> s == sym, pauli_ops) - 1
symboltoint(i::Integer) = i
inttosymbol(int::Integer) = pauli_ops[int+1]
inttosymbol(s::Symbol) = s

function symboltoint(oper::Vector{Symbol})
    nq = length(oper)
    intoper = getinttype(nq)(0)
    for (ii, symb) in enumerate(oper)
        intoper = setelement!(intoper, ii, symboltoint(symb))
    end
    return intoper
end

function inttosymbol(int::Integer, n_qubits::Integer)
    symbs = [:I for _ in 1:n_qubits]
    for ii in 1:n_qubits
        symbs[ii] = inttosymbol(getelement(int, ii))
    end
    return symbs
end


import Base.show
function show(op::Unsigned)
    max_qubits_in_integer = round(Int, bitsize(typeof(op)) / 2)
    show(op, max_qubits_in_integer)
end

function show(op::Unsigned, n::Int)
    max_qubits_in_integer = round(Int, bitsize(typeof(op)) / 2)
    nind = min(max_qubits_in_integer, n)
    print_string = prod("$(inttosymbol(getelement(op, ii)))" for ii in 1:nind)
    println(print_string)

end

function show(d::Dict; max_qubits_shown=16)
    ## TODO
    println("Not yet implemented.")
end
