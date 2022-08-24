
@testset "Basics" begin
    # Ensure basic functionality:
    # - things multiply correctly
    # - things add correctly
    # - things throw errors at appropiate times

    using SparseArrays, LinearAlgebra

    @testset "Equality" begin
        @test XX() == XX(1)
        @test YY(2) == 2.0*YY()
        @test ZZ(4) == 2.0*ZZ(2)
        @test X() != Y()
        @test Z([1,2]) == [1,2]*Z()
        @test 2*Hopp() != 3*Hopp()
        @test [1,0]*Y() != [0,1]*Y()
        @test Z([1]) != Z(1)
        @test_throws ErrorException [1,2] * ZZ() # 1-body information for 2-body term
        @test_throws ErrorException rand(3,4) * ZZ() # non-square interaction matrix
    end

    @testset "Add" begin
        @testset "Terms" begin
            @test XX(1)+YY(2) isa SpinModels.SumOfTerms
            @test X([2,3]) + Y([1,2]) isa SpinModels.SumOfTerms
            @test_throws ErrorException X([1,2]) + Y(3) # one has spatial structure, the other doesn't
            @test_throws ErrorException XX([1,2]) + YY([1,2,3]) # different spatial structure
            @test [0 1; 1 0] * XX() + [1,2]*Z() isa SpinModels.SumOfTerms
            @test_throws ErrorException [0 1; 1 0] * XX() + [1,2,3]*Z() # different #spins
        end

        @testset "SumOfTerms" begin
            @test (XX(1) + YY(2)) + ZZ(3) isa SpinModels.SumOfTerms
            @test_throws ErrorException (X(1) + Y(2)) + Z([3])
            @test [4]*(X(1) + Y(2)) + Z([3]) isa SpinModels.SumOfTerms
            @test [4]*(X(1) + Y(2)) + Z([3]) isa SpinModels.SumOfTerms

            @test (XX(1) + YY(2)) + (X(3) + Y(4)) isa SpinModels.SumOfTerms # note: Maybe this should throw? You can never multiply it with something to then convert to a matrix, so it seems pretty useless...
            @test_throws ErrorException (XX(1) + YY(2)) + [3]*(X(3) + Y(4))
            @test_throws ErrorException [0 1; 1 0]*(XX(1) + YY(2)) + (X(3) + Y(4))
            @test_throws ErrorException [0 1; 1 0]*(XX(1) + YY(2)) + [3]*(X(3) + Y(4))
            @test [0 1; 1 0]*(XX(1) + YY(2)) + [3,4]*(X(3) + Y(4)) isa SpinModels.SumOfTerms
        end
    end

    @testset "Multiply" begin
        @test X(1) == X()*1
        @test 2*X(1) == X(2)
        @test 2*X(2) == X(4)
        @test [1,2]*X(2) == 2*X([1,2])
        @test (X()*[1,2])*3 == X([3,6])
        @test_throws ErrorException [1,2]*X([1,2])
    end

    @testset "Matrix conversion" begin

        𝟙 = I(2)
        σx = [0 1; 1 0]
        σy = [0 -im; im 0]
        σz = Diagonal([1,-1])

        𝟙𝟙 = I(4)
        xx = kron(σx,σx)
        yy = kron(σy,σy)
        zz = kron(σz,σz)
        hopp = sparse([2,3], [3,2], [1,1], 4,4)

        @testset "Small matrices" begin
            @test sparse([1]*X()) == σx
            @test sparse([1]*Y()) == σy
            @test sparse([1]*Z()) == σz

            J = [0 1; 0 0]
            @test sparse(J * XX()) == xx
            @test sparse(J * YY()) == yy
            @test sparse(J * ZZ()) == zz
            @test sparse(J * Hopp()) == hopp

            # Wrapper types
            @test sparse(J' * XX()) == xx
            @test sparse(Matrix(J') * XX()) == xx
            @test sparse((J'+J)/2 * XX()) == xx
            @test sparse(sparse(J) * XX()) == xx
            @test sparse(Symmetric(J) * XX()) == 2*xx
            @test sparse(Hermitian(J) * XX()) == 2*xx
            @test sparse(UpperTriangular(J) * XX()) == xx
            @test sparse(Tridiagonal(J) * XX()) == xx
            @test sparse(SymTridiagonal(Symmetric(J)) * XX()) == 2*xx
        end

        @testset "slightly bigger matrices" begin
            h = [1,2.0]
            for (op, σ) in zip((X,Y,Z), (σx,σy,σz))
                @test sparse(h*op()) ≈ h[1]*kron(σ, 𝟙) + h[2]*kron(𝟙, σ) # todo change if basis ordering changes!
                @test sparse(h*op()) ≈ sparse([h[1], 0]*op() + [0, h[2]]*op())
                @test sparse(op(h)) ≈ sparse(h*op())
            end

            J1 = 1; J2= 2.0;
            J = sparse([1,3], [2,4],[J1,J2], 4,4) # 2 uncoupled pairs of spins
            for (op, σ) in zip((XX,YY,ZZ,Hopp), (xx,yy,zz,hopp)) # XXX
                @test sparse(J*op()) ≈ J1*kron(σ, 𝟙𝟙) + J2*kron(𝟙𝟙, σ) # todo change if basis ordering changes!
                @test sparse(op(J)) ≈ sparse(J*op())
            end
        end

        @testset "Composition" begin
            h = rand(5)
            J = rand(5,5)
            @test sparse(J*XX() + h*Z()) ≈ sparse(J*XX()) + sparse(h*Z())
            @test sparse(2*(J*XX() + h*Z())) ≈ 2*sparse(J*XX()) + 2*sparse(h*Z())
            @test sparse(J*XX() + J*ZZ()) ≈ sparse(J*XX()) + sparse(J*ZZ())
            @test sparse(h*X() + h*Z()) ≈ sparse(h*X()) + sparse(h*Z())
        end
    end
end
