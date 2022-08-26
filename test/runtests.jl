using Test
using SpinModels

@testset "SpinModels.jl" begin
    include("basics.jl")
    include("geometry.jl")
    include("interactions.jl")
    include("xxznumerics-consistency.jl")

    @testset "NN compat" begin
        @test interaction_matrix(PowerLaw(5), NN(Chain(5))) ≈ interaction_matrix(NN(PowerLaw(5)), Chain(5))
        @test interaction_matrix(PowerLaw(5), NN(Chain(5))) ≈ interaction_matrix(NN(PowerLaw(5)), NN(Chain(5)))
    end
end
