# SpinModels.jl

**Construct all of your favorite spin ``\frac{1}{2}`` models with ease.**

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
TODO

## Full documentation
```@contents
Pages = ["full-docs/hamiltonian.md","full-docs/geometry.md","full-docs/interaction.md"]
Depth = 1
```