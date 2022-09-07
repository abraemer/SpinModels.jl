abstract type Hamiltonian end

struct SumOfTerms <: Hamiltonian
    termlist#::Vector{Term}
end

# make kind a type parameter instead?
# Then could change SumOfTerms to hold a tuple and be generic as well
# no more type instability! and loops over terms could probably be inlined by the compiler...
# Let's benchmark that at some point :)
struct Term{T} <: Hamiltonian
    kind::Symbol
    coefficient::T
end

const ONEBODY_TERMS = [:X, :Y, :Z]
const TWOBODY_TERMS = [:XX, :YY, :ZZ, :Hopp, :XXX]

# Convenience methods
for (sym, sym_name) in ((:X, "\\sigma_x^{(i)}"),
                        (:Y, "\\sigma_y^{(i)}"),
                        (:Z, "\\sigma_z^{(i)}"))
    @eval begin
        @doc """
            $($sym)(J=1)

        Represents a ``$($sym_name)`` term. Multiply with a vector h to get the term
        ``\\sum_i h_i $($sym_name) .``
        """
        $sym(J=1) = Term($(Meta.quot(sym)), J)
    end
end

for (sym, sym_name) in ((:XX, "\\sigma_x^{(i)} \\sigma_x^{(j)}"),
                        (:YY, "\\sigma_y^{(i)} \\sigma_y^{(j)}"),
                        (:ZZ, "\\sigma_z^{(i)} \\sigma_z^{(j)}"),
                        (:Hopp, "\\sigma_+^{(i)} \\sigma_-^{(j)} + \\mathrm{h.c.}"))
    @eval begin
        @doc """
            $($sym)(J=1)

        Represents a ``$($sym_name)`` term. Multiply with a matrix J to get the term
        ``\\sum_{i\\neq j} J_{ij} $($sym_name) .``

        Note: Two-body terms do 'double-counting' of the couplings meaning the sum runs over `i,j` independently and not `i < j`.
        """
        $sym(J=1) = Term($(Meta.quot(sym)), J)
    end
end

"""
    FlipFlop(J=1)

Same as [`Hopp`](@ref).
"""
FlipFlop = Hopp

"""
    XXZ(Δ, J=1)

Represents the terms: `Hopp(2*J)+Δ*ZZ(J)`
"""
XXZ(Δ, J=1) = Hopp(2*J) + ZZ(Δ*J)

Base.:(==)(t1::Term, t2::Term) = t1.kind == t2.kind && t1.coefficient == t2.coefficient
Base.hash(t::Term, h) = hash(t.coefficient, hash(t.kind, h))
# TODO comparison of SumOfTerms
# complication: order of terms should not matter
# implement after deciding how to solve this
# possibly after putting the termtypes into SumOfTerm type and enforcing some order among terms

"""
    _nbodyness(term)

Number of spins the term acts on.
"""
function _nbodyness end
_nbodyness(term) = _nbodyness(Val(term.kind))
for sym in ONEBODY_TERMS
    @eval _nbodyness(::$(typeof(Val(sym)))) = 1
end
for sym in TWOBODY_TERMS
    @eval _nbodyness(::$(typeof(Val(sym)))) = 2
end

"""
    _check_array(dims, nbody)

Bad name. Raise an error if an array with size `dims` is not suitable as coefficient for an `nbody` term.
Generally a n-body term needs a n-dimensional array which axes are all equal.
"""
function _check_array(dims, nbody)
	length(dims) == nbody || error("Incompatible coefficient array with size $(dims) for $nbody-body term!")
	all(==(dims[1]), dims) || error("Coefficient array with size $(dims) has unequal dimensions!")
	true
end

# multiplication/division of Term with array/number -> just affects the coefficient field
# Note that an array determines the spatial structure, i.e. which spins the operator acts on.
# As such we can only multiply by an array a single time!
Base.:*(n::Number, t::Term) = Term(t.kind, n*t.coefficient)
Base.:*(t::Term, n::Number) = Term(t.kind, t.coefficient*n)
Base.:/(t::Term, n::Number) = Term(t.kind, t.coefficient / n)

Base.:*(a::AbstractArray, t::Term{<:Number}) = _check_array(size(a), _nbodyness(t)) && Term(t.kind, a*t.coefficient)
Base.:*(t::Term{<:Number}, a::AbstractArray) = a*t

Base.:*(::AbstractArray, ::Term{<:AbstractArray}) = error("Term already has a spatial structure! You cannot multiply an array again.")
Base.:*(::Term{<:AbstractArray}, ::AbstractArray) = error("Term already has a spatial structure! You cannot multiply an array again.")

## Adding Terms
# Adding 2 terms gives a SumOfTerms, which just contains both.
# TODO: Merge terms of same kind??
# TODO: Consistent ordering of terms? -> When using Terms as type param this would reduce #compilations
"""
   _check_add_compatible(term1, term2)

Two terms can be added, if either:
1. both have *no* spatial information
2. both carry spatial information for the same number of spins
Otherwise error.
"""
function _check_add_compatible end
_check_add_compatible(t1::Term{<:AbstractArray}, t2::Term{<:AbstractArray}) = nspins(t1) == nspins(t2) || error("Couplings shapes give different #spins: $(size(t1.coefficient)) vs $(size(t2.coefficient))")
_check_add_compatible(t1::Term{<:Number}, t2::Term{<:Number}) = true
_check_add_compatible(t1::Term, t2::Term) = error("Terms cannot be added due to inconsistent spatial information! $(size(t1.coefficient)) vs. $(size(t2.coefficient))")

Base.:+(t1::Term, t2::Term) = _check_add_compatible(t1,t2) && SumOfTerms([t1,t2])
Base.:+(t1::Term, t2::SumOfTerms) = _check_add_compatible(t1,t2.termlist[1]) && SumOfTerms([t1; t2.termlist])
Base.:+(t1::SumOfTerms, t2::Term) = _check_add_compatible(t2,t1.termlist[1]) && SumOfTerms([t1.termlist; t2])
Base.:+(t1::SumOfTerms, t2::SumOfTerms) = _check_add_compatible(t1.termlist[1],t2.termlist[1]) && SumOfTerms([t1.termlist; t2.termlist])

# simple subtraction
Base.:-(t::Hamiltonian) = -1 * t
Base.:-(t1::Hamiltonian, t2::Hamiltonian) = t1 + (-t2)

# multiplication/division acts on each term separately
Base.:*(x::Union{Number, AbstractArray}, t::SumOfTerms) = SumOfTerms(Ref(x) .* t.termlist)
Base.:*(t::SumOfTerms, x) = SumOfTerms(t.termlist .* Ref(x))
Base.:/(t::SumOfTerms, x) = SumOfTerms(t.termlist ./ Ref(x))

function nspins end
"""
    nspins(term)
    nspins(sum_of_terms)

Return the number of spins in the underlying geometry if that information was already provided.
"""
nspins(::Hamiltonian) = error("Not implemented! This should never happen!")
nspins(::Term{<:Number}) = error("No spatial structure defined yet!")
nspins(t::Term{<:AbstractArray}) = size(t.coefficient, 1)
nspins(sot::SumOfTerms) = nspins(sot.termlist[1])

## Conversion to Matrix
## Right now this counts spins from the LEFT meaning J[L] <-> 1st spins and J[1]<->Lth spin
# I probably want to change this...
# Can I abstract that to a macro and let the user choose which ordering he wants?

LinearAlgebra.Matrix(H::Hamiltonian) = Matrix(SparseArrays.sparse(H))
SparseArrays.SparseMatrixCSC(H::Hamiltonian) = SparseArrays.sparse(H)

# Convenience method to avoid code duplication
function _termlist end
_termlist(t::Term) = [t]
_termlist(sot::SumOfTerms) = sot.termlist

function SparseArrays.sparse(H::Hamiltonian)
    terms = _termlist(H)
    Ls = _coo_length_estimate.(terms)
    L = sum(Ls)
    I = zeros(Int,L)
    J = zeros(Int,L)
    V = zeros(_eltype_estimate(H), L)

    offset = 0
    for (L, term) in zip(Ls, terms)
        @views construct_COO!(I[1+offset:L+offset],J[1+offset:L+offset],V[1+offset:L+offset], term)
        offset += L
    end
    SparseArrays.sparse(I,J,V, 2^nspins(H), 2^nspins(H))
end

_eltype_estimate(t::Term) = t.kind == :Y ? complex(eltype(t.coefficient)) : eltype(t.coefficient)
_eltype_estimate(sot::SumOfTerms) = promote_type(_eltype_estimate.(sot.termlist)...)

_coo_length_estimate(term::Term) = _coo_length_estimate(Val(term.kind), term)
_coo_length_estimate(::Union{Val{:ZZ}, Val{:Z}}, term) = 2^nspins(term) # diagonal terms
_coo_length_estimate(::Any, term) = _count_nonzeros(term.coefficient)*2^nspins(term)
_coo_length_estimate(::Val{:Hopp}, term) = _count_nonzeros(term.coefficient) * 2^(nspins(term)-1)
_coo_length_estimate(::Val{:XXX}, term) = (_count_nonzeros(term.coefficient)+2) * 2^(nspins(term)-1)

_count_nonzeros(a) = binomial(size(a,1), length(size(a)))
##
## More precise estimate only make sense when the COO-generating code also respects the zeros!
##
#_count_nonzeros(a) = min(binomial(size(a,1), length(size(a))), _count_nonzeros_specialized(a))
# _count_nonzeros_specialized(a::Any) = length(a)
# _count_nonzeros_specialized(a::SparseArrays.AbstractSparseVector) = SparseArrays.nnz(a)
# _count_nonzeros_specialized(a::SparseArrays.AbstractSparseMatrix) = SparseArrays.nnz(a) # TODO Symmetry?
# # unwrap LinearAlgebra's wrapper types
# _count_nonzeros_specialized(a::Union{LinearAlgebra.Symmetric, LinearAlgebra.Hermitian, LinearAlgebra.UpperTriangular, LinearAlgebra.LowerTriangular, LinearAlgebra.Transpose}) = _count_nonzeros_specialized(parent(a))
# _count_nonzeros_specialized(a::LinearAlgebra.Diagonal) = 0 # Diagonal values don't count...
# _count_nonzeros_specialized(a::Union{LinearAlgebra.Bidiagonal,LinearAlgebra.Tridiagonal,LinearAlgebra.SymTridiagonal}) = size(a,1)-1 # off-diagonal entries

construct_COO!(I,J,V, t::Term) = _construct_COO!(I,J,V, Val(t.kind), t.coefficient)

# ZZ
function _construct_COO!(I,J,V, ::Val{:ZZ}, coeff)
    N = size(coeff,1)
    I .= 1:2^N
    J .= 1:2^N
    V .= 0
    for i in 1:N
        for j in i+1:N
            _zcorrelator_values!(V, N, coeff[j,i]+coeff[i,j], i, j)
        end
    end
end

"""
	    _zcorrelator_values!(V, N, Jij, i, j)

*ADD* to `V` the diagonal of Jᵢⱼ σzⁱσzʲ.

# Parameters
- `V` array of length 2^`N`. Will be mutated!
- `N` number of sites in the system
- `Jij` interaction strength
- `i`,`j` sites to couple
"""
function _zcorrelator_values!(V, N, Jij, i, j)
    # 𝟙(sec1) ⊗ σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    #             | 𝟙(sec2) |     0    |   | 𝟙(sec3) |    0     |
    # = 𝟙(sec1) ⊗ | ------- | -------- | ⊗ | ------- | -------- |
    #             |    0    | -𝟙(sec2) |   |    0    | -𝟙(sec3) |


    #             |           [ 𝟙(sec3) |     0    ] |                                  |
    #             | 𝟙(sec2) ⊗ [ ------- | -------- ] |              0                   |
    #             |           [    0    | -𝟙(sec3) ] |                                  |
    # = 𝟙(sec1) ⊗ | -------------------------------- | -------------------------------- |
    #             |                                  |           [ -𝟙(sec3) |    0    ] |
    #             |              0                   | 𝟙(sec2) ⊗ [ -------- | ------- ] |
    #             |                                  |           [    0     | 𝟙(sec3) ] |
    sec1 = 2^(i-1)
    sec2 = 2^(j-i-1)
    sec3 = 2^(N-j)

    blocksize3 = 2sec3 # σz ⊗ 𝟙(sec3)
    blocksize23 = 2sec2*blocksize3 # size of: σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    #sec1 loop
    for offset1 in blocksize23 .* (0:sec1-1)
        # upper block
        for offset2 in offset1 .+ blocksize3 .* (0:sec2-1)
            V[1 .+ offset2 .+ (0:sec3-1)] .+= Jij
            V[1 .+ offset2 .+ (sec3:2sec3-1)] .+= -Jij
        end
        # lower block
        for offset2 in offset1 .+ blocksize3 .* (sec2:2sec2-1)
            V[1 .+ offset2 .+ (0:sec3-1)] .+= -Jij
            V[1 .+ offset2 .+ (sec3:2sec3-1)] .+= Jij
        end
    end
end

## Hopping/FlipFlop

function _construct_COO!(A,B,V, ::Val{:Hopp}, J)
    N = size(J,1)
    fill!(V, 1)
    length = 2^(N-1)
    at = 0
    for i in 1:N
        for j in i+1:N
            range = at .+ (1:length)
            @views _hopping_coordinates!(A[range], B[range], N, i, j)
            V[range] .*= J[j, i]+J[i, j]
            at += length
        end
    end
    ## _hopping values! does the transpose
    ## no measureable performance difference
    #copyto!(A,at+1,B,1,at)
    #copyto!(B,at+1,A,1,at)
    #copyto!(V,at+1,V,1,at)
end

"""
    _hopping_coordinates!(A, B, N, i, j)

Compute the coordinates of the matrix entries of the hopping operator between site i and J
in a spin chain. The values are 1 everywhere, so there is no input/output for these in this
function.

Hopping operator: σ₊ⁱσ₋ʲ + σ₋ⁱσ₊ʲ = 𝟙⊗…⊗𝟙⊗σ₊⊗𝟙⊗…⊗𝟙⊗σ₋⊗𝟙⊗…⊗𝟙 + transpose

# Parameters
 - `A`,`B` coordinate vectors. These should have length 2^(N-2) will be overwritten!
 - `N` number of sites in the system
 - `i`,`j` sites the excitation can hop between. 1 ≤ i,j ≤ N
"""
function _hopping_coordinates!(A,B,N,i,j)
    i > j && return _hopping_coordinates!(A,B,N,j,i)
    # 𝟙(sec1) ⊗ σ+ ⊗ 𝟙(sec2) ⊗ σ- ⊗ 𝟙(sec3) + transpose
    #             |                                  |                                  |
    #             |                0                 |                 0                |
    #             |                                  |                                  |
    # = 𝟙(sec1) ⊗ | -------------------------------- | -------------------------------- |
    #             |           [    0    |  𝟙(sec3) ] |                                  |
    #             | 𝟙(sec2) ⊗ [ ------- | -------- ] |                 0                |
    #             |           [    0    |    0     ] |                                  |
    sec1 = 2^(i-1)
    sec2 = 2^(j-i-1)
    sec3 = 2^(N-j)

    blocksize3 = 2sec3 # σz ⊗ 𝟙(sec3)
    blocksize23 = 2sec2*blocksize3 # size of: σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    currentIndex = 0
    #sec1 loop
    for offset1 in blocksize23 .* (0:sec1-1)
        for offset2 in blocksize3 .* (0:sec2-1)
            rowOffset2 = offset1 + offset2 # left columns # offset for the row-COORDINATE
            colOffset2 = offset1 + offset2 + blocksize3*sec2 # bottom rows # offset for the col-COORDINATE

            A[currentIndex .+ (1:sec3)] .= colOffset2 .+ (1:sec3) # upper rows
            B[currentIndex .+ (1:sec3)] .= rowOffset2 .+ (sec3+1:2sec3) # right columns
            currentIndex += sec3

            ## now transpose
            B[currentIndex .+ (1:sec3)] .= colOffset2 .+ (1:sec3) # upper rows
            A[currentIndex .+ (1:sec3)] .= rowOffset2 .+ (sec3+1:2sec3) # right columns
            currentIndex += sec3
        end
    end
end

## Z

function _construct_COO!(I,J,V, ::Val{:Z}, coeff)
    N = size(coeff,1)
    I .= 1:2^N
    J .= 1:2^N
    V .= 0
    for i in 1:N
        _z_field!(V, N, coeff[i], i)
    end
end

"""
    _z_field!(V, N, hz, i)

*ADD* to `V` the diagonal of hz σzⁱ.

# Parameters
 - `V` array of length 2^`N`. Will be mutated!
 - `N` number of sites in the system
 - `hz` field strength
 - `i` site to act on
"""
function _z_field!(V, N, hz, i)
    # σz^(i) = 𝟙(sec1) ⊗ σz ⊗  𝟙(sec2)
    sec1 = 2^(i-1)
    sec2 = 2^(N-i)

    blocksize2 = 2sec2 # σz ⊗ 𝟙(sec2)

    for offset1 in blocksize2 .* (0:sec1-1)
        V[offset1 .+ (1:sec2)] .+= hz
        V[offset1 .+ (sec2+1:2sec2)] .-= hz
    end
end

## X and Y
function _construct_COO!(I,J,V, ::Val{:X}, coeff)
    N = size(coeff,1)
    length = 2^N
    at = 0
    for i in 1:N
        range = at .+ (1:length)
        V[range] .= coeff[i]
        @views _xy_coordinates!(I[range], J[range], N, i)
        at += length
    end
end

function _construct_COO!(I,J,V, ::Val{:Y}, coeff)
    # Y = im * [0 -1; 1 0]
    N = size(coeff,1)
    length = 2^N
    at = 0
    for i in 1:N
        range = at .+ (1:length)
        @views _xy_coordinates!(I[range], J[range], N, i)
        for ind in range
            input_spin = 1 & ((J[ind]-1) >> (N-i)) # 0 ↔ spin up, 1↔spin down
            V[ind] = coeff[i] * im * (1 - 2*(input_spin)) # +im if input state is up ≡ bit not set
        end
        at += length
    end
end

function _xy_coordinates!(I, J, N, i)
    # σx^(i) = 𝟙(sec1) ⊗ σx ⊗  𝟙(sec2)
    I .= xor.(0:(2^N)-1, 2^(N-i)) .+ 1 # flip ith spin in column
    J .= 1:(2^N)
end

## XX YY

function _construct_COO!(I,J,V, ::Val{:XX}, coeff)
    N = size(coeff,1)
    # I[:] .= 1:2^N
    # J[:] .= 1:2^N
    length = 2^N
    at = 0
    for i in 1:N
        for j in i+1:N
            range = at .+ (1:length)
            @views _xxyy_coordinates!(I[range], J[range], N, i, j)
            V[range] .= coeff[j, i]+coeff[i, j]
            at += length
        end
    end
end

function _construct_COO!(I,J,V, ::Val{:YY}, coeff)
    N = size(coeff,1)
    length = 2^N
    at = 0
    for i in 1:N
        for j in i+1:N
            range = at .+ (1:length)
            @views _xxyy_coordinates!(I[range], J[range], N, i, j)
            for ind in range
                input_state = J[ind]-1
                input_spin1 = 1 & (input_state >> (N-i)) # 0 ↔ spin up, 1↔spin down # N-i -> i-1 for other ordering!
                input_spin2 = 1 & (input_state >> (N-j)) # 0 ↔ spin up, 1↔spin down
                V[ind] = (coeff[j,i]+coeff[i,j]) * (2*xor(input_spin1,input_spin2) - 1) # +1 if input spins are anti-aligned, -1 if aligned
            end
            at += length
        end
    end
end

function _xxyy_coordinates!(I, J, N, i,j)
    # σx^(i) = 𝟙(sec1) ⊗ σx ⊗  𝟙(sec2)
    I .= xor.(xor.(0:(2^N)-1, 2^(N-i)), 2^(N-j)) .+ 1 # flip ith and jth spins in column
    J .= 1:(2^N)
end

## XXX
function _construct_COO!(I,J,V, ::Val{:XXX}, coeff)
    lengthZZ = 2^size(coeff,1)
    @views _construct_COO!(I[1:lengthZZ],J[1:lengthZZ],V[1:lengthZZ], Val(:ZZ), coeff)
    @views _construct_COO!(I[lengthZZ:end],J[lengthZZ:end],V[lengthZZ:end], Val(:XX), coeff)
end
