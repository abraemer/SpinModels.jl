### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 956589f2-13f2-11ed-230e-55379de914d0
# ╠═╡ show_logs = false
using PlutoUI,SparseArrays,LinearAlgebra

# ╔═╡ 1bc4f3fa-c5ef-44ec-8d27-696595db2e9a
html"""<style>main { max-width: 40%;}</style>"""

# ╔═╡ 3ab30833-694d-4490-8ff1-2059f4208480
begin
	abstract type Hamiltonian end

	struct SumOfTerms <: Hamiltonian
	    termlist#::Vector{Term}
	end
	
	struct Term{T} <: Hamiltonian
		kind::Symbol
	    coefficient::T
	end
	
	struct OneBodyTerm{T} <: Hamiltonian
	    coefficient
	end
	

	for sym in (:XX,:YY,:ZZ,:HOPP,:X,:Y,:Z)
		@eval $sym(J=1) = Term($(Meta.quot(sym)), J)
	end
end

# ╔═╡ 518c0c80-0ee3-4731-a470-acfa2df93b8d
begin
	function _nbodyness end
	_nbodyness(t::Term) = _nbodyness(Val(t.kind))
	for sym in (:X,:Y,:Z)
		@eval _nbodyness(::$(typeof(Val(sym)))) = 1
	end
	for sym in (:XX,:YY,:ZZ,:HOPP)
		@eval _nbodyness(::$(typeof(Val(sym)))) = 2
	end
end

# ╔═╡ 823028d8-e0b7-45d9-9f75-ee5ecf79940e
function _check_array(dims, nbody)
	length(dims) == nbody || error("Incompatible coefficient array with size $(dims) for $nbody-body term!")
	all(==(dims[1]), dims) || error("Coefficient array with size $(dims) has unequal dimensions!")
	true
end

# ╔═╡ a1bdd4e7-141e-45b4-bb54-f3d32217f217
begin
	Base.:*(n::Number, t::Term) = Term(t.kind, n*t.coefficient)
	Base.:*(t::Term, n::Number) = Term(t.kind, t.coefficient*n)
	Base.:*(a::AbstractArray, t::Term{<: Number}) = _check_array(size(a), _nbodyness(t)) && Term(t.kind, a*t.coefficient)
	Base.:*(t::Term{<: Number}, a::AbstractArray) = a*t
	Base.:/(t::Term, n::Number) = Term(t.kind, t.coefficient / n)
end

# ╔═╡ a6d45ac2-c66a-46cb-a41a-d2b2c6f2cddc
XX()/2

# ╔═╡ 00eb02bd-552c-4471-85d0-2799aee82913
begin
	nspins(t::Term{<:Number}) = error("No spatial structure defined yet!")
	nspins(t::Term{<:AbstractArray}) = size(t.coefficient, 1)
	nspins(sot::SumOfTerms) = nspins(sot.termlist[1])
end

# ╔═╡ 3ce9de97-c3e3-443c-98f6-626730247fb3
begin
	_add_compatible(t1::Term, t2::Term) = false
	_add_compatible(t1::Term{<:AbstractArray}, t2::Term{<:AbstractArray}) = nspins(t1) == nspins(t2) || error("Couplings shapes give different #spins: $(size(t1.coefficient)) vs $(size(t2.coefficient))")
	_add_compatible(t1::Term{<:Number}, t2::Term{<:Number}) = true
	_check_add_compatible(t1,t2) = _add_compatible(t1,t2) || error("Terms cannot be added due to inconsistent spatial information! $(size(t1.coefficient)) vs. $(size(t2.coefficient))")
	
	Base.:+(t1::Term, t2::Term) = _check_add_compatible(t1,t2) && SumOfTerms([t1,t2])
	Base.:+(t1::Term, t2::SumOfTerms) = 
		_check_add_compatible(t1,t2.termlist[1]) && SumOfTerms([t1; t2.termlist])
	Base.:+(t1::SumOfTerms, t2::Term) =
		_check_add_compatible(t2,t1.termlist[1]) && SumOfTerms([t1.termlist; t2])
	Base.:+(t1::SumOfTerms, t2::SumOfTerms) =
		_check_add_compatible(t1.termlist,t2.termlist[1]) && SumOfTerms([t1.termlist; t2.termlist])
	Base.:*(x, t::SumOfTerms) = SumOfTerms(Ref(x) .* t.termlist)
	Base.:*(t::SumOfTerms, x) = SumOfTerms(t.termlist .* Ref(x))
end

# ╔═╡ c346d83a-12cc-4427-80e1-c773c53d1f5d
2*Term(:XX, 1)

# ╔═╡ 81dc12fe-44dc-4470-9159-2765d10efe55
2*(rand(10,10)*XX())

# ╔═╡ b0ffa57a-ce4e-4298-90a5-619705c66fcf
begin
	Base.:-(t::Hamiltonian) = -1 * t
	Base.:-(t1::Hamiltonian, t2::Hamiltonian) = t1 + (-t2)
end

# ╔═╡ 6596aed9-6a7f-46b7-b6e8-d37fb408c02b
XX() - YY()

# ╔═╡ 05240c9f-6b78-443d-a0ef-0c7c6e375f16
-1 .* ([XX(),YY()])

# ╔═╡ 40e4534f-5ea1-46bf-94ed-fe91744bedf4
rand(4,4) * (XX()+YY() + -0.7*ZZ())

# ╔═╡ 0c555bc2-f856-495d-a2e9-924ec7bbcf67
rand(4,4) * (XX()+YY()) + -0.7*ZZ()

# ╔═╡ c07d4b92-5ed7-4ccd-b559-c3b59d7aeda1
rand(4,4) * (XX()+YY()) + rand(5,5)*-0.7*ZZ()

# ╔═╡ cc5f6047-b6b6-4d35-a78c-e1ba6aff857b
rand(5,4) * HOPP()

# ╔═╡ 6ceafa48-ec77-4621-8ae2-6f971f0f77fc
typeof([2.0*XX(),YY()])

# ╔═╡ efebd035-948c-44aa-9597-715ce21c47a5
supertypes(SparseVector)

# ╔═╡ 8ed1d7ef-7c00-406d-acdd-d1b9393fcd9a
supertypes(SparseMatrixCSC)

# ╔═╡ 16f71b33-ff5f-465f-b8f5-059f10a248dd
_eltype_estimate(::Hamiltonian) = Float64 #  TODO Y Term

# ╔═╡ 1009d6de-e795-4c77-b979-0876a5a9a93b
begin
	_coo_length_estimate(term::Term) = _coo_length_estimate(Val(term.kind), term)

	# diagonal terms
	_coo_length_estimate(::Union{Val{:ZZ}, Val{:Z}}, term) = 2^nspins(term)
	# other terms
	_coo_length_estimate(::Val{:HOPP}, term) = 
	let N = _count_nonzeros(term.coefficient)
		N*2^(nspins(term)-1)
	end
	_coo_length_estimate(::Any, term) = 
	let N = _count_nonzeros(term.coefficient)
		N*2^nspins(term)
	end

	function _count_nonzeros end
	_count_nonzeros(a::DenseVector) = length(a)
	_count_nonzeros(a::AbstractMatrix) = binomial(size(a,1), length(size(a))) # Todo: higher dimensions?
	_count_nonzeros(a::AbstractSparseVector) = nnz(a)
	_count_nonzeros(a::AbstractSparseMatrix) = nnz(a)÷2 # TODO Symmetry?
end

# ╔═╡ 52db2ed3-e9d7-4b2e-b27f-b34518c9246a
begin
	"""
	    _construct_COO!(A,B,V, ::Val{:HOPP}, J)
	
	Generate the COO vectors for ∑ᵢⱼ Jᵢⱼ (σ₊ⁱσ₋ʲ + σ₋ⁱσ₊ʲ)
	
	# Parameters
	 - `A`,`B`,`V` vectors for column (`A`) and row indices (`B`) respectively values (`V`)
	 - `J` interaction matrix of size `N`x`N`. Assumed symmetric!
	
	 If `J` is `N`x`N`, then `A`,`B` and `V` need to have length (N-1)N*2^(N-2)
	"""
	function _construct_COO!(A,B,V, ::Val{:HOPP}, J)
		N = size(J,1)
	    fill!(V, 1)
	    length = 2^(N-1)
	    at = 0
	    for i in 1:N
	        for j in i+1:N
	            range = at .+ (1:length)
	            @views _hopping_coordinates!(A[range], B[range], N, i, j)
	            V[range] .*= J[j, i]
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
end

# ╔═╡ e3ebe672-fdf3-4fc6-b3bb-0aef2ca5d2e8
begin
	"""
    _construct_COO!(I, J, V, ::Val{:ZZ}, C)

	Generate the COO vectors of ∑ᵢⱼCᵢⱼ σzⁱσzʲ.
	
	# Parameters
	 - `I`,`J`,`V` vectors for column (`A`) and row indices (`B`) respectively values (`V`)
	 - `C` interaction matrix of size `N`x`N`
	
	`I`,`J` and `V` need to have length 2^N.
	"""
	function _construct_COO!(I,J,V, ::Val{:Z}, coeff)
		N = size(coeff,1)
	    I[:] .= 1:2^N
	    J[:] .= 1:2^N
	    V[:] .= 0
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
end

# ╔═╡ 41085b33-f48f-4509-bc98-d0209acbfb6e
begin
	construct_COO!(I,J,V, t::Term) = _construct_COO!(I,J,V, Val(t.kind), t.coefficient)
	
	"""
    _construct_COO!(I, J, V, ::Val{:ZZ}, C)

	Generate the COO vectors of ∑ᵢⱼCᵢⱼ σzⁱσzʲ.
	
	# Parameters
	 - `I`,`J`,`V` vectors for column (`A`) and row indices (`B`) respectively values (`V`)
	 - `C` interaction matrix of size `N`x`N`
	
	`I`,`J` and `V` need to have length 2^N.
	"""
	function _construct_COO!(I,J,V, ::Val{:ZZ}, coeff)
		N = size(coeff,1)
	    I[:] .= 1:2^N
	    J[:] .= 1:2^N
	    V[:] .= 0
	    for i in 1:N
	        for j in i+1:N
	            _zcorrelator_values!(V, N, coeff[j,i], i, j)
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
end

# ╔═╡ 3861b6ec-c019-4f09-91d3-b214d6b5d714
begin
	function SparseArrays.sparse(term::Term)
		L = _coo_length_estimate(term)
		I = zeros(Int,L)
		J = zeros(Int,L)
		V = zeros(_eltype_estimate(term), L)
		construct_COO!(I,J,V, term)
		sparse(I,J,V, 2^nspins(term), 2^nspins(term))
	end

	function SparseArrays.sparse(sot::SumOfTerms)
		Ls = _coo_length_estimate.(sot.termlist)
		L = sum(Ls)
		I = zeros(Int,L)
		J = zeros(Int,L)
		V = zeros(_eltype_estimate(sot), L)

		offset = 0
		for (L, term) in zip(Ls, sot.termlist)
			@views construct_COO!(I[1+offset:L+offset],J[1+offset:L+offset],V[1+offset:L+offset], term)
			offset += L
		end
		sparse(I,J,V, 2^nspins(sot), 2^nspins(sot))
	end
	
	LinearAlgebra.Matrix(H::Hamiltonian) = Matrix(sparse(H))
end

# ╔═╡ ed618a22-8e08-48b6-ab42-ecf129df5105
let J = [0 1;1 0],
	H = J*(HOPP() + -0.5*ZZ()) + [-0.05, 0.05] * Z()
	@show H
	sparse(H)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
PlutoUI = "~0.7.39"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═1bc4f3fa-c5ef-44ec-8d27-696595db2e9a
# ╠═956589f2-13f2-11ed-230e-55379de914d0
# ╠═3ab30833-694d-4490-8ff1-2059f4208480
# ╠═518c0c80-0ee3-4731-a470-acfa2df93b8d
# ╠═823028d8-e0b7-45d9-9f75-ee5ecf79940e
# ╠═a1bdd4e7-141e-45b4-bb54-f3d32217f217
# ╠═c346d83a-12cc-4427-80e1-c773c53d1f5d
# ╠═81dc12fe-44dc-4470-9159-2765d10efe55
# ╠═3ce9de97-c3e3-443c-98f6-626730247fb3
# ╠═b0ffa57a-ce4e-4298-90a5-619705c66fcf
# ╠═05240c9f-6b78-443d-a0ef-0c7c6e375f16
# ╠═40e4534f-5ea1-46bf-94ed-fe91744bedf4
# ╠═0c555bc2-f856-495d-a2e9-924ec7bbcf67
# ╠═c07d4b92-5ed7-4ccd-b559-c3b59d7aeda1
# ╠═cc5f6047-b6b6-4d35-a78c-e1ba6aff857b
# ╠═a6d45ac2-c66a-46cb-a41a-d2b2c6f2cddc
# ╠═6596aed9-6a7f-46b7-b6e8-d37fb408c02b
# ╠═6ceafa48-ec77-4621-8ae2-6f971f0f77fc
# ╠═00eb02bd-552c-4471-85d0-2799aee82913
# ╠═3861b6ec-c019-4f09-91d3-b214d6b5d714
# ╠═efebd035-948c-44aa-9597-715ce21c47a5
# ╠═8ed1d7ef-7c00-406d-acdd-d1b9393fcd9a
# ╠═16f71b33-ff5f-465f-b8f5-059f10a248dd
# ╠═1009d6de-e795-4c77-b979-0876a5a9a93b
# ╠═41085b33-f48f-4509-bc98-d0209acbfb6e
# ╠═52db2ed3-e9d7-4b2e-b27f-b34518c9246a
# ╠═ed618a22-8e08-48b6-ab42-ecf129df5105
# ╠═e3ebe672-fdf3-4fc6-b3bb-0aef2ca5d2e8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
