## Test the datatypes module with all possible constructors and adders

using Test

#TODO: Add example tests for PauliString constructors
function createpaulistring(nq)
    symbol = rand([:I, :X, :Y, :Z])
    qind = rand(1:nq)
    coeff = randn()
    PauliString(nq, symbol, qind, coeff)

    symbols = rand([:I, :X, :Y, :Z], min(nq, 4))
    qinds = rand(1:nq, min(nq, 4))
    coeff = randn()
    pstr = PauliString(nq, symbols, qinds, coeff)
    print(pstr)
    return pstr
end

function createpaulisum(nq)
    PauliSum(nq)

    pstr = createpaulistring(nq)
    PauliSum(nq, pstr)

    pstr = createpaulistring(nq)
    psum = PauliSum(pstr)
    print(psum)
    return psum
end

function addtopaulisum(nq)
    psum = createpaulisum(nq)
    pstr = createpaulistring(nq)
    pstr_temp = psum + pstr
    @test isa(pstr_temp, PauliSum)
    pstr_temp = pstr + pstr
    @test isa(pstr_temp, PauliSum)
    add!(psum, pstr)
    @test getcoeff(psum, pstr.term) == pstr.coeff

    symbol = rand([:I, :X, :Y, :Z])
    qind = rand(1:nq)
    coeff = randn()
    add!(psum, symbol, qind, coeff)
    @test getcoeff(psum, symbol, qind) == coeff

    symbols = rand([:I, :X, :Y, :Z], min(nq, 4))
    qinds = rand(1:nq, min(nq, 4))
    coeff = randn()
    psum2 = createpaulisum(nq)
    add!(psum2, symbols, qinds, coeff)
    @test getcoeff(psum2, symbols, qinds) == coeff

    psum3 = psum + psum2
    print(topaulistrings(psum3))
    psum2 - psum3

    return psum
end

@testset "PauliString Tests" begin
    pstr = createpaulistring(7)
    wrapped_pstr = wrapcoefficients(pstr, PauliFreqTracker)
end

# Test PauliSum from Dict creation
function test_paulisum_from_dict()
    psum = PauliSum(3, Dict([:I, :I, :I] => 1.5, [:I, :I, :Y] => 1.0))

    # Collect and return keys and values for testing
    return collect(keys(psum.terms)), collect(values(psum.terms))
end

# Test subtraction of PauliSum
function subtractpaulisums()
    psum1 = PauliSum(3, Dict([:I, :I, :I] => 1.5, [:I, :I, :Y] => 1.0))
    psum2 = PauliSum(3, Dict([:I, :I, :I] => 1.5))

    return psum1 - psum2
end

# Test overloading methods for PauliSum
# TODO(YT): Add tests for *, /, + overloading for PauliSum

@testset "PauliSum Tests" begin
    # Subtest for PauliSum from Dict
    paulis, pauli_cs = test_paulisum_from_dict()
    @test paulis == [symboltoint([:I, :I, :Y]), symboltoint([:I, :I, :I])]
    @test pauli_cs == [1.0, 1.5]

    # Subtest for subtracting PauliSum
    result_psum = subtractpaulisums()
    expected_psum = PauliSum(3, Dict([:I, :I, :Y] => 1.0))
    @test result_psum == expected_psum
end
