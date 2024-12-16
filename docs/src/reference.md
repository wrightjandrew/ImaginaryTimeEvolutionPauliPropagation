# Reference

```@contents
Pages = ["reference.md"]
```

## Gates

```@docs
Gate
ParametrizedGate
StaticGate

PauliRotation

CliffordGate
clifford_map
createcliffordmap
reset_clifford_map!

TGate

ParametrizedNoiseChannel
PauliNoise
DepolarizingNoise
DephasingNoise
AmplitudeDampingNoise
PauliXNoise
PauliYNoise
PauliZNoise

FrozenGate
freeze
```

## Circuits

```@docs
countparameters
bricklayertopology
staircasetopology
get2dtopology
get2dstaircasetopology
ibmeagletopology
hardwareefficientcircuit
efficientsu2circuit
tfitrottercircuit
tiltedtfitrottercircuit
heisenbergtrottercircuit
su4ansatz
qcnnansatz
appendSU4!
```

## Pauli Types

```@docs
PauliString
PauliSum
PauliStringType
PauliType
getpauli
setpauli
symboltoint
PauliPropagation.add!
PauliPropagation.set!
```

## Pauli Algebra

```@docs
commutes
commutator
pauliprod
countweight
countxy
countyz
containsXorY
containsYorZ
```

## PathProperties

```@docs
PathProperties
NumericPathProperties
numcoefftype
tonumber
wrapcoefficients
```

## Propagate

```@docs
propagate
propagate!
applymergetruncate!
applytoall!
applyandadd!
mergeandempty!
```

## State Overlap

```@docs
overlapbyorthogonality
overlapwithzero
overlapwithplus
overlapwithcomputational
overlapwithmaxmixed
overlapwithpaulisum
PauliPropagation.filter
PauliPropagation.filter!
zerofilter
zerofilter!
plusfilter
plusfilter!
```

## Apply

```@docs
apply
```

## Surrogate

```@docs
NodePathProperties
evaluate!
reset!
```


## Index
```@index
Pages = ["reference.md]
```