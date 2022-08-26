@testset "interactions" begin

    using Random, Statistics, LinearAlgebra
    rng = Random.Xoshiro(5) # my favorite seed

    N = 5
    chain = Chain(N)
    dists = distance_matrix(chain)

    constInt = ConstantInteraction(3)
    @test constInt(chain) == interaction_matrix(constInt, chain)
    @test constInt(chain) == 3*(dists .> 0)

    powerlaw2 = PowerLaw(2)
    ref = 1 ./ distance_matrix(chain) .^ 2
    ref[diagind(ref)] .= 0
    @test interaction_matrix(powerlaw2, chain) ≈ ref
    @test interaction_matrix(powerlaw2, chain) .^ 2 ≈ interaction_matrix(PowerLaw(4), chain)

    @testset "NNInteraction $geom with $int" for geom in (Chain(5), PBC(Chain(5))),
                                                  int in (PowerLaw(5), ConstantInteraction(3))
        nn = NN(int)
        @test nn isa SpinModels.NNInteraction
        @test_throws MethodError SpinModels._compute_interaction(nn, geom, [1], [2]) # not defined

        full_int = int(geom)
        nn_mask = full_int .>= 1
        @test nn(geom) == nearest_neighbour_from_interactions(full_int)
        @test nn(geom)[nn_mask] == full_int[nn_mask]
        @test all(nn(geom)[.!nn_mask] .== 0)
    end
end
