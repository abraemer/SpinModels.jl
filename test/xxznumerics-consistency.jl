#
# Check generated matrices by comparing them to xxznumerics
#
@testset "XXZNumerics.jl consistency" begin
    import Pkg
    Pkg.add(url="https://github.com/abraemer/XXZNumerics.jl")
    import XXZNumerics
    using LinearAlgebra, Random, SparseArrays
    Ns = [7,8,9,10,11]
    Random.seed!(5)
    σx = XXZNumerics.σx
    σy = XXZNumerics.σy
    σz = XXZNumerics.σz

    checkprecision(x) = @test sum(abs2, x) < 1e-15

    @testset "1-body terms" begin
        @testset "$opname" for (opname, xxznum_op, hm_op) in (("σx", σx, X()), ("σy", σy, Y()),("σz", σz, Z()))
            for N in Ns, i in 1:10 # 10 repetitions
                J = rand(N)
                checkprecision(sum(J[i]*XXZNumerics.single_spin_op(xxznum_op, i, N) for i in 1:N) .- sparse(J*hm_op))
            end
        end
    end

    @testset "2-body terms" begin

        @testset "$opname" for (opname, xxznum_op, hm_op) in (("σxσx", σx, XX()), ("σyσy", σy, YY()),("σzσz", σz, ZZ()))
            for N in Ns, i in 1:10 # 10 repetitions
                J = rand(N,N)
                J = J'J # symmetrize
                J[diagind(J)] .= 0

                checkprecision(sum([2*J[i,j]*XXZNumerics.correlator(xxznum_op, i, j, N) for i in 1:N for j in i+1:N]) .- sparse(J*hm_op))
            end
        end

        @testset "hopping" begin
            for N in Ns, i in 1:10 # 10 repetitions
                J = rand(N,N)
                J = J'J # symmetrize
                J[diagind(J)] .= 0

                checkprecision(sum([J[i,j]*(XXZNumerics.correlator(σx, i, j, N)+XXZNumerics.correlator(σy, i,j, N)) for i in 1:N for j in i+1:N]) .- sparse(J*Hopp()))
            end
        end
    end

    @testset "xxzmodel" begin
        for N in Ns, i in 1:10 # 10 repetitions
            J = rand(N,N)
            J = J'J # symmetrize
            J[diagind(J)] .= 0
            Δ = randn()

            checkprecision(XXZNumerics.xxzmodel(J, Δ) .- sparse(J/8*XXZ(Δ))) # /4 to convert pauli->spin-op, /2 due to double counting
        end
    end

    # more tests for compat between hopping, xxzmodel and so on
    # Note to self: Currently Hopp/FlipFlop = XX() + YY() but in XXZNumerics theres a factor 1/2!
end
