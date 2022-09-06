# Coupling matrices

The task of creating a coupling matrix is divided into two parts:
1. The [Geometry](@ref) describes how many spins in what spatial configuration there are
2. The [Interaction](@ref) knows how to translate the positions into coupling strengths

So generally you want to first build up the desired geometry and then apply the correct interaction to that. Of course you can always directly provide a coupling matrix if these standards do not cover your needs. Just keep in mind the [index order](@ref index_order).

## Geometry
There are several simple base geometries: [`Chain`](@ref), [`NoisyChain`](@ref) and [`Box`](@ref). The latter two don't admit to fixed positions and instead draw them randomly each time [`positions`](@ref) is called on them.

Then there are several modifiers, you can apply:
 - [`PBC`](@ref) to enforce periodic boundary conditions. This basically changes how distances are computed within the geometry.
 - [`Blockaded`](@ref) to enforce a minimal distance between spins. This only makes sense to apply to disordered geometries.
 - [`NN`](@ref) to only keep nearest neighbour distances.

You may use function like [`positions`](@ref), [`distance_matrix`](@ref) and [`nearest_neighbor_from_distances`](@ref) on geometry objects.

## Interaction
The next step is to define the desired interaction type. Currently implemented are [`ConstantInteraction`](@ref) and [`PowerLaw`](@ref). You can also modify these to be nearest-neighbor only with [`NN`](@ref).

Once you have defined your desired geometry and interaction you can get the final interaction matrix via [`interaction_matrix`](@ref) (or equivalently by calling the interaction object with the geometry).