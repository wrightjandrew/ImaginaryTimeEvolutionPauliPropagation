### miscgates.jl
##
# A file for gates that don't fit into the other categories.
##
###



struct TGate <: StaticGate
    qind::Int

    @doc """
        TGate(qind::Integer)

    Returns a T gate acting on qubit `qind`.
    It acts on qubit `qind` like a `PauliRotation(:Z, qind)` with angle π/4.
    """
    function TGate(qind::Integer)
        if qind <= 0
            throw(ArgumentError("Qubit index must be a positive integer."))
        end
        return new(qind)
    end
end