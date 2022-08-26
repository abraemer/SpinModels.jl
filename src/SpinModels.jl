module SpinModels
using Random
import LinearAlgebra
import SparseArrays

export X,Y,Z,XX,YY,ZZ,Hopp,FlipFlop,XXX,XXZ
include("basics.jl")

export nspins, positions, distance_matrix, nearest_neighbour_from_distances, nearest_neighbour_from_distances!
export Box, Chain, NoisyChain, Blockaded, PBC
include("geometry.jl")

export interaction_matrix, nearest_neighbour_from_interactions, nearest_neighbour_from_interactions!
export ConstantInteraction, PowerLaw
include("interactions.jl")

export NN
NN(geom::Geometry, k=1) = NNGeometry(geom; k)
NN(int::Interaction, k=1) = NNInteraction(int; k)

end #module
