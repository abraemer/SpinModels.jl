# Quick reference

## Terms
Add terms together to define the structure of your Hamiltonian and multiply with constants and coupling matrices (only once!) to set up a concrete realization.
- [`X`](@ref)`(h)`, [`Y`](@ref)`(h)` and [`Z`](@ref)`(h)` stand for ``\sum_i h_i \sigma_\alpha^{(i)},`` where ``\alpha`` is ``x``, ``y`` or ``z``.
- [`XX`](@ref)`(J)`, [`YY`](@ref)`(J)` and [`ZZ`](@ref)`(J)` stand for ``\sum_{i,j} J_{i,j} \sigma_\alpha^{(i)}\sigma_\alpha^{(j)},`` where ``\alpha`` is ``x``, ``y`` or ``z``.
- [`FlipFlop`](@ref)`(J)` (or equivalently [`Hopp`](@ref)`(J)`) is the same as `0.5*J*(XX()+YY())`
- [`XXZ`](@ref)`(Δ, J)` is the same as `2*Hopp(J)+Δ*ZZ(J)` which is the same as `XX(J)+YY(J)+Δ*ZZ(J)`

## Geometries
Current implemented are (`N` always denotes the total number of spins):
- [`Chain`](@ref)`(N)`
- [`Box`](@ref)`(N, dims)`
- [`NoisyChain`](@ref)`(N, σ; spacing=1)`
- [`PartiallyFilledChain`](@ref)`(numspins, numsites; spacing=1)`
- [`Blockaded`](@ref)`(geometry; retries=1000, blockade=1.0)`
- [`PBC`](@ref)`(geometry)`
- [`NN`](@ref)`(geometry, k=1)`

!!! note

    The order of [`Blockaded`](@ref), [`PBC`](@ref) and [`NN`](@ref) is irrelevant. They sort them themselves.

## Interactions
Currently implemented are
- [`ConstantInteraction`](@ref)`(value)`
- [`PowerLaw`](@ref)`(α)`
- [`NN`](@ref)`(interaction, k=1)`

!!! note

    There are QoL overloads on Geometry that apply `ConstantInteraction(1)` automatically. Thus you should probably not need to use it directly. Please give feedback if that does unexpected things.

## Note about nearest neighbor (`NN`)
For isotropic interactions (i.e. all currently implemented ones) it does not matter whether you apply `NN` to the interaction or the geometry. In principle for anisotropic interaction there will be a subtle difference:
- `NN` on geometry is based on distances and will remove all but the `k`smallest distances
- `NN` on interactions is based on coupling strength and will remove all but the `k`strongest couplings