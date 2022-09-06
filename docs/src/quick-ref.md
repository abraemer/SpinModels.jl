# Quick reference

## Terms
- [`X`](@ref)`(h)`, [`Y`](@ref)`(h)` and [`Z`](@ref)`(h)` stand for ``\sum_i h_i \sigma_\alpha^{(i)},`` where ``\alpha`` is ``x``, ``y`` or ``z``.
- [`XX`](@ref)`(J)`, [`YY`](@ref)`(J)` and [`ZZ`](@ref)`(J)` stand for ``\sum_{i,j} J_{i,j} \sigma_\alpha^{(i)}\sigma_\alpha^{(j)},``
where ``\alpha`` is ``x``, ``y`` or ``z``.
- [`FlipFlop`](@ref)`(J)` (or equivalently [`Hopp`](@ref)`(J)`) is the same as `0.5*J*(XX()+YY())`
- [`XXZ`](@ref)`(Δ, J)` is the same as `2*Hopp(J)+Δ*ZZ(J)` which is the same as `XX(J)+YY(J)+Δ*ZZ(J)`

## Geometries
Current implemented are (`N` always denotes the total number of spins):
- [`Chain`](@ref)`(N)`
- [`Box`](@ref)`(N, dims)`
- [`NoisyChain`](@ref)`(N, σ, spacing=1)`
- [`Blockaded`](@ref)`(geometry; retries=1000, blockade=1.0)`
- [`PBC`](@ref)`(geometry)`
- [`NN`](@ref)`(geometry, k=1)`

## Interactions
Currently implemented are
- [`ConstantInteraction`](@ref)`(value)`
- [`PowerLaw`](@ref)`(α)`
- [`NN`](@ref)`(interaction, k=1)`

## Note about nearest neighbor (`NN`)
