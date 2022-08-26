### Helpers for distance computation
function _euclidean(p1, p2)
    return sqrt(sum(@. (p1 - p2)^2 ))
end

function _euclidean_pbc(p1, p2, lengths)
    diffs = _euclidean_pbc.(p1, p2, lengths)
    return sqrt(sum(diffs.^2))
end

function _euclidean_pbc(x1::Number, x2::Number, length::Number)
    diff = abs(mod(x1, length) - mod(x2, length))
    return diff > length/2 ? length-diff : diff
end
############################################
################# Geometry #################
############################################

abstract type Geometry end

### Public interface
function nspins end
function positions end
function distance_matrix end

### Private interface
function _distance end

### Default implementations
distance_matrix(geom; kwargs...) = distance_matrix(geom, positions(geom; kwargs...))

function distance_matrix(geom, positions)
    # consistency check?
    # @assert nspins(geom) == size(positions, 2)
    N = size(positions,2)
    res = zeros(N, N)
    for i in 1:N
        for j in i+1:N
            dist = _distance(geom, positions[:,i], positions[:,j])
            res[i,j] = dist
            res[j,i] = dist
        end
    end
    return res
end

############################################
#### BaseGeometry

abstract type BaseGeometry <: Geometry end

### Implementation: Chain
struct Chain{T} <: BaseGeometry
    L::Int
    spacing::T
    Chain(L; spacing=1) = new{typeof(spacing)}(L,spacing)
end

nspins(c::Chain) = c.L
positions(c::Chain) = (1:nspins(c))' .* c.spacing
_distance(::Chain, p1, p2) = _euclidean(p1,p2)


### Implementation: NoisyChain
struct NoisyChain{T} <: BaseGeometry
    L::Int
    spacing::T
    σ::Float64
    NoisyChain(L, σ; spacing=1) = new{typeof(spacing)}(L,spacing, σ)
end

nspins(c::NoisyChain) = c.L
positions(c::NoisyChain; rng=Random.GLOBAL_RNG) = ((1:nspins(c)) .* c.spacing .+ c.σ .* (rand(rng, nspins(c)) .- 0.5))'
_distance(::NoisyChain, p1, p2) = _euclidean(p1,p2)


### Implementation: Box
struct Box{T} <: BaseGeometry
    N::Int
    dims::Vector{T}
end

nspins(c::Box) = c.N
positions(b::Box; rng=Random.GLOBAL_RNG) = rand(rng, length(b.dims), nspins(b)) .* b.dims
_distance(::Box, p1, p2) = _euclidean(p1,p2)

############################################
#### GeometryModifier

abstract type GeometryModifier{G<:Geometry} <: Geometry end

### Private Interface
function _parent end

### Default implementation
nspins(gm::GeometryModifier) = nspins(_parent(gm))
positions(gm::GeometryModifier; kwargs...) = positions(_parent(gm); kwargs...)
_distance(hm::GeometryModifier, p1, p2) = _distance(_parent(hm), p1, p2)
_parent(gm::GeometryModifier) = gm.geometry


### Implementation: Blockaded
# retries < 0 -> try until success
struct Blockaded{G} <: GeometryModifier{G}
    geometry::G
    retries::Int
    blockade::Float64
    Blockaded(geometry; retries=1000, blockade=1.0) = new{typeof(geometry)}(geometry, retries, blockade)
end

## helpers
_lower_diag_indices(mat) = (CartesianIndex(i,j) for j in 1:size(mat,2) for i in j+1:size(mat,1))
_each_lower_diag_entry(mat) = (mat[I] for I in _lower_diag_indices(mat))
# assumes symmetry of distances
_fulfills_blockade_condition(distances, blockade) = all(>(blockade), _each_lower_diag_entry(distances))

## Overwrite defaults
function distance_matrix(blockaded::Blockaded, positions)
    distances = distance_matrix(blockaded.geometry, positions)
    if !_fulfills_blockade_condition(distances, blockaded.blockade)
        minI = argmin(I->distances[I], _lower_diag_indices(distances))
        error("""Given positions do not satisfy the blockade condition. \
                 Distance at $minI is $(distances[minI]) \
                 but blockade radius is $(blockaded.blockade)""")
    end
    return distances
end

function positions(block::Blockaded; kwargs...)
    # retries < 0 -> try until success
    tries = block.retries
    while tries != 0
        tries -= 1
        pos = positions(block.geometry; kwargs...)
        distances = distance_matrix(block.geometry, pos)
        if _fulfills_blockade_condition(distances, block.blockade)
            return pos
        end
    end
    error("Blockade condition could not be fulfilled in $(block.retries) tries.")
end


### Implementation: NNGeometry
struct NNGeometry{G<:Geometry} <: GeometryModifier{G}
    geometry::G
    k::Int
    NNGeometry(geometry; k=1) = new{typeof(geometry)}(geometry, k)
end

## Helpers
nearest_neighbour_from_distances(distance_mat, k=1) = nearest_neighbour_from_distances!(copy(distance_mat), k)

function nearest_neighbour_from_distances!(distance_mat, k=1)
    # select k+1 smallest distances
    # k+1 because spin to itself always has dist 0
    d_cutoff = [sort!(unique(col))[k+1] for col in eachcol(distance_mat)]
    distance_mat[distance_mat .> d_cutoff'] .= Inf # not nearest nearest_neighbour = infinitely far away
    # could also consider converting to sparse but probably not worth it most of the times?
    # for sparse matrices simply remove the entry and not set it to Inf?
    return distance_mat
end

## Overwrite defaults
distance_matrix(nng::NNGeometry, positions) = nearest_neighbour_from_distances!(distance_matrix(_parent(nng), positions), nng.k)


### Implementation: PBC
struct PBC{G<:BaseGeometry} <: GeometryModifier{G}
    geometry::G
end

# switch PBC to the inside for consistency
PBC(b::Blockaded) = Blockaded(PBC(_parent(b)); retries=b.retries, blockade=b.blockade)
PBC(nng::NNGeometry) = NNGeometry(PBC(_parent(nng)); k=nng.k)

## Overwrite defaults
_distance(pbcchain::PBC{<:Union{Chain, NoisyChain}}, p1, p2) = _euclidean_pbc(p1[1], p2[1], nspins(pbcchain.geometry)*pbcchain.geometry.spacing)
_distance(boxpbc::PBC{<:Box}, p1, p2) = _euclidean_pbc(p1, p2, boxpbc.geometry.dims)
