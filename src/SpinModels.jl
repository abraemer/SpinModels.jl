module SpinModels
using Random
import LinearAlgebra
import SparseArrays

export X,Y,Z,XX,YY,ZZ,Hopp,FlipFlop,XXX,XXZ
export nspins
include("basics.jl")

export positions, distance_matrix, nearest_neighbor_from_distances, nearest_neighbor_from_distances!
export Box, Chain, NoisyChain, Blockaded, PBC
include("geometry.jl")

export interaction_matrix, nearest_neighbor_from_interactions, nearest_neighbor_from_interactions!
export ConstantInteraction, PowerLaw
include("interactions.jl")

export NN

"""
    NN(geometry, k=1)

Modify the geometry to only keep the `k` closest connections from each spin. Please see
[`nearest_neighbor_from_distances!`](@ref) for further information.

!!! note

    This is **not** identical to [`NN(::SpinModels.Interaction)`](@ref) if the interaction is anisotropic.
"""
NN(geom::Geometry, k=1) = NNGeometry(geom; k)

"""
    NN(interaction, k=1)

Modify the interaction to only keep the `k` strongest couplings from each spin. Please see
[`nearest_neighbor_from_interactions!`](@ref) for further information.

!!! note

    This is **not** identical to [`NN(::SpinModels.Geometry)`](@ref) if the interaction is anisotropic.
"""
NN(int::Interaction, k=1) = NNInteraction(int; k)

end #module
