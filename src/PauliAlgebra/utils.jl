
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


inttostring(op::Unsigned) = prod("$(inttosymbol(getelement(op, ii)))" for ii in 1:Int(bitsize(op) / 2))

import Base.show
function show(op::Integer)
    println(inttostring(op))

end

function show(op::Integer, n::Int)
    max_qubits_in_integer = round(Int, bitsize(typeof(op)) / 2)
    nind = min(max_qubits_in_integer, n)

    print_string = inttostring(op)[1:nind]
    println(print_string)

end


function show(d::Dict; max_lines=20)
    show(d, Int(bitsize(first(d)[1]) / 2); max_lines=max_lines)

end

function show(d::Dict, n::Int; max_lines=20)
    header = "$(typeof(d)) with $(length(d)) entries:"
    println(header)
    for (ii, (op, coeff)) in enumerate(d)
        if ii > max_lines
            println("  ⋮")
            break
        end
        println("  ", inttostring(op), " => ", coeff)
    end

end