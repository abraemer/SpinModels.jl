# SpinModels.jl

**Construct all of your favorite spin-``\frac{1}{2}`` models with ease.**

## Recipe
1. Define a coupling matrix like `J = PowerLaw(6)(Chain(10))`
2. Add the relevant terms of your Hamiltonian and apply couplings like `H = J*(XX() + YY() + Δ*ZZ())`
3. Convert to sparse/dense with `SparseArrays.sparse(H)`/`LinearAlgebra.Matrix(H)`

## Usage examples
The (standard) MBL Hamiltonian ``\hat{H} = J\sum_i \left(\sigma_+^{(i)}\sigma_-^{(j)} + \sigma_-^{(i)}\sigma_+^{(j)} + \Delta \sigma_z^{(i)}\sigma_z^{(j)}\right) + \sum_i h_i \sigma_z^{(i)}`` could be constructed like so:
```julia
N = 10 # number of spins
h = 5.0 # strength of disorder
Δ = 1.0 # strength of Ising coupling, here XXX model
geometry = NN(Chain(N)) # nearest neighbor chain...
J_matrix = interaction_matrix(ConstantInteraction(J), geometry) # ... with constant interactions
# could also use ConstantInteraction(J)(geometry)
h_vec = h*(rand(N) .- 1) # random fields
H = J_matrix * (Hopping() + 2*Δ*ZZ()) + h*Z()
```

This object ```H``` can now be converted to either sparse or dense matrices using `SparseArrays.sparse(H)` or `LinearAlgebra.Matrix(H)` respectively.

## Structure of this package
A spin model's Hamiltonian consists of the sum of different terms, such as `XX` or `Z`, that denote different Pauli operator combinations (here: ``\sigma^{(i)}_x \sigma^{(j)}_x`` and ``\sigma_z^{(i)}``). To be able to realize the Hamiltonian's matrix, each term needs to have information about the coupling strengths. For one-body terms this is a vector and for two-body terms a matrix. The sizes of these coupling arrays correspond to the amount of spins in the model (and as such need to be equal across all terms).

For ease of use, there are several helpful structs and functions implemented to make constructing commonly used coupling types easy. They are organized into [Geometry](@ref)s and [Interaction](@ref)s. See [Coupling matrices](@ref).

## [A word on ordering the basis states](@id index_order)
The basis of a spin-``\frac{1}{2}`` model can be easily represented by numbers in binary, where each ``0`` denotes a ``|\uparrow\rangle`` and each ``1`` denotes a ``|\downarrow\rangle``. The 3-spin state ``|\uparrow\downarrow\uparrow\rangle`` for example gets the number ``010_2 = 2_{10}`` and as such is the third basis state (1-based counting).

There are now two ways of organizing the basis with respect to the coupling matrices:
The coupling ``J_{i,j}`` couples the ``i``th spin with the ``j``th spin, where the ``i``th spin is
1. the spin at position ``i`` in the state vector (order like ``|1\rangle\otimes|2\rangle\otimes|3\rangle\ldots``)
2. the spin corresponding to the ``i``th digit in the binary expansion with value ``2^{i-1}`` (order like ``\ldots|3\rangle\otimes|2\rangle\otimes|1\rangle``)

I'd argue that for the sake of doing bit manipulation on the state numbers the latter is more convenient as one does not need to know the total amount of spins. Currently this library counts its spins from the **LEFT** corresponding to the **FORMER** (1.) option (due to legacy code reasons).

Ultimately this issue boils down to the inconsistency in writing generally from left to right except when writing numbers, where the lowest valued place is on the right and increasing to the left.

## Full documentation
```@contents
Pages = ["full-docs/hamiltonian.md","full-docs/geometry.md","full-docs/interaction.md"]
Depth = 1
```