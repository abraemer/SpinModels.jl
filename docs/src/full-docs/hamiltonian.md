# Full docs for Hamiltonian
## Types
These `Term` types may be added together freely. Multiplication with scalars works as expected. You can only impose spatial information (via multiplication with a vector or matrix respective of the type of term) once! You can only add together compatible terms in terms of spatial structure.
```@docs
X
Y
Z
XX
YY
ZZ
Hopp
FlipFlop
XXZ
```


## Functions
```@docs
nspins(::SpinModels.Hamiltonian)
```