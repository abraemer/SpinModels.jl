using SpinModels
using Test

@testset "SpinModels.jl" begin
    include("basics.jl")
    include("xxznumerics-consistency.jl")
end
