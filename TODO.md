# SpinModels

## Todo implementation
 - decide what to do about basis ordering
 - compare performance to more type-stable version
 - Add automated benchmarks (for learning purposes)?
 - Think about automatic meanfield functions
 - Register package?
 - ╔ (Optimization) Matrix construction
 - ║ Generalize matrix construction to directly construct dense matrices, build up sparse matrices
 - ║ Include possibility to directly construct Matrix in symmetrized subspace -> SpinSymmetry.jl
 - ║ Code of building: setindex!(arr,x,y,v)
 - ║ For dense -> set directly, for sparse -> have a COO and just append to I,J,V vectors
 - ║ For symmetry sector: put x,y through the projection matrix and apply factor to v, then setindex!
 - ╚ 

## Todo Documentation
 - ~~Explain basis ordering! (after deciding)~~
 - MORE Usage + examples ~~(after implementing the geometry stuff)~~
 - benchmarks vs. kron?
 - ~~Mention SpinSymmetry.jl~~

