@testset "geometry" begin

    using Random, Statistics, LinearAlgebra
    rng = Random.Xoshiro(5) # my favorite seed

    @testset "Chain" begin
        c1 = Chain(5)
        c2 = Chain(5;spacing=2)
        reference_positions = reshape([1,2,3,4,5], 1, 5)

        @test nspins(c1) == 5
        @test nspins(c2) == 5

        @test positions(c1) ≈ reference_positions
        @test positions(c2) ≈ 2*positions(c1)

        @test distance_matrix(c1) ≈ abs.(reference_positions .- reference_positions')
        @test distance_matrix(c2) ≈ 2*distance_matrix(c1)
    end

    @testset "NoisyChain" begin
        # compare no noise to normal chain
        nc1 = NoisyChain(5, 0)
        nc2 = NoisyChain(5, 0; spacing=2)
        c1 = Chain(5)
        c2 = Chain(5;spacing=2)

        @test nspins(nc1) == 5
        @test nspins(nc2) == 5

        @test positions(nc1) ≈ positions(c1)
        @test positions(nc2) ≈ positions(c2)

        @test distance_matrix(nc1) ≈ distance_matrix(c1)
        @test distance_matrix(nc2) ≈ distance_matrix(c2)


        σ = 2
        nc3 = NoisyChain(5, σ)
        reference_positions = positions(c1)
        for i in 1:20
            pos = positions(nc3; rng)
            relpos = pos .- reference_positions
            @test all(-σ/2 .< relpos) && all(relpos .< σ/2)
        end

        # check if average positions are also close to no noise chains
        # variance of a single position should be 2/3/σ (σ/2)^2
        meanpos = mean(positions(nc3; rng) for i in 1:10000)
        varpos = mean(abs2.(positions(nc3; rng) .- reference_positions) for i in 1:10000)
        σexpected = 2/(3σ)*(σ/2)^2
        # @info mean(meanpos .- reference_positions)
        # @info mean(varpos)
        @test mean(meanpos .- reference_positions) < σ/2 / 10000 / nspins(nc3)
        @test 0.9*σexpected < mean(varpos) < 1.1*σexpected
    end

    @testset "PartiallyFilledChain" begin
        for i in -4:4
            @test_throws ErrorException PartiallyFilledChain(5,i)
        end
        @test_throws ErrorException PartiallyFilledChain(-5,5)

        # compare no noise to normal chain
        pfc1 = PartiallyFilledChain(5, 5)
        pfc2 = PartiallyFilledChain(5, 5; spacing=2)
        c1 = Chain(5)
        c2 = Chain(5;spacing=2)

        @test nspins(pfc1) == 5
        @test nspins(pfc2) == 5

        @test positions(pfc1) ≈ positions(c1)
        @test positions(pfc2) ≈ positions(c2)

        @test distance_matrix(pfc1) ≈ distance_matrix(c1)
        @test distance_matrix(pfc2) ≈ distance_matrix(c2)

        for (N,L) in [(10,10),(10,50),(10,200),(10,1000)]
            geom = PartiallyFilledChain(N,L)
            @test all(allunique, positions(geom) for _ in 1:100)
        end

        SHOTS = 1000
        pfc3 = PartiallyFilledChain(5, 10)
        data = sort!(mapreduce(vec, vcat, positions(pfc3) for _ in 1:SHOTS))
        @test 1:10 == unique(data) # all positions where generated
        counts = diff(unique(i->data[i], eachindex(data)))
        reference_value = SHOTS*pfc3.numspins/pfc3.numsites
        deviation = sqrt(reference_value) # shot noise is ~√N
        @test all(reference_value-2deviation .< counts .< reference_value+2deviation)
    end

    @testset "Box" begin
        box1 = Box(2, [1])
        box2 = Box(3, [1,2.0])
        box3 = Box(4, [1,2,3])

        @test nspins(box1) == 2
        @test nspins(box2) == 3
        @test nspins(box3) == 4

        @test size(positions(box1; rng)) == (1,2)
        @test size(positions(box2; rng)) == (2,3)
        @test size(positions(box3; rng)) == (3,4)

        shots = 10000
        meanpos1 = mean(mean(positions(box1; rng); dims=2) for i in 1:shots)
        meanpos2 = mean(mean(positions(box2; rng); dims=2) for i in 1:shots)
        meanpos3 = mean(mean(positions(box3; rng); dims=2) for i in 1:shots)

        # @info meanpos1 .- [0.5]
        # @info meanpos2 .- [0.5,1]
        # @info meanpos3 .- [0.5,1,1.5]

        # safety factor of 3 on the expected deviations
        @test sum(abs, meanpos1 .- [0.5]) < 3 * 0.5/√shots
        @test sum(abs, meanpos2 .- [0.5,1]) < 3 * norm([0.5,1])/√shots
        @test sum(abs, meanpos3 .- [0.5,1,1.5]) < 3 * norm([0.5,1,1.5])/√shots
    end

    @testset "Blockaded" begin
        @test PBC(Blockaded(Chain(5))) isa Blockaded # normalizes Blockaded to the outside

        N = 17
        chain = Chain(N; spacing = 2)
        b_chain = Blockaded(chain)
        @test positions(chain) == positions(b_chain)
        @test distance_matrix(chain) == distance_matrix(b_chain)

        @test_throws ErrorException positions(Blockaded(Chain(5; spacing = 0.5))) # can't generate blockaded positions
        @test_throws ErrorException distance_matrix(Blockaded(Chain(2)), [1 1.5]) # positions are already blockaded
    end

    @testset "NNGeometry" begin
        for nn in (NN(Chain(5)), NN(PBC(Chain(5))))
            @test nn isa SpinModels.NNGeometry
            @test positions(nn) == positions(nn.geometry)
            full_dists = distance_matrix(nn.geometry)
            nnmask = full_dists .<= 1
            @test distance_matrix(nn)[nnmask] == full_dists[nnmask]
            @test all(distance_matrix(nn)[.!nnmask] .== Inf)
        end
    end

    @testset "PBC" begin
        for N in (16,17) # even/odd
            c = Chain(N)
            pbcc = PBC(c)
            @test nspins(pbcc) == nspins(c)
            @test positions(pbcc) == positions(c)

            pbc_dists = distance_matrix(pbcc)
            @test allequal(circshift(col, -i+1) for (i,col) in enumerate(eachcol(pbc_dists))) # translation invariance
            @test all(pbc_dists[i,1]==pbc_dists[N-i+2,1] for i in 2:div(N,2)) # spatial reflection

            pbc_nc = PBC(NoisyChain(N, 0))
            @test distance_matrix(pbcc) ≈ distance_matrix(pbc_nc)

            pbc_box = PBC(Box(2,[1]))
            shots = 10000
            mean_dist = mean(distance_matrix(pbc_box; rng)[1,2] for i in 1:shots)
            @test mean_dist ≈ 0.25 atol=3/√shots # safety factor of 3
        end
        # propagation checks
        @test PBC(Chain(5)) isa PBC
        @test PBC(NN(Chain(5))).geometry isa PBC
        @test PBC(Blockaded(NN(Chain(5)))).geometry.geometry isa PBC
        @test PBC(NN(Blockaded(Chain(5)))).geometry.geometry isa PBC
    end
end
