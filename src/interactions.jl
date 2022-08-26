###########################################
############### Interaction ###############
###########################################

abstract type Interaction end

### Private Interface
function _compute_interaction end

### Public Interface
function interaction_matrix end

### default implementation
(interaction::Interaction)(args...;kwargs...) = interaction_matrix(interaction, args...; kwargs...)

interaction_matrix(interaction, geom; kwargs...) = interaction_matrix(interaction, geom, positions(geom; kwargs...))

function interaction_matrix(interaction, geom, positions)
    N = nspins(geom)
    res = zeros(N, N)
    dists = distance_matrix(geom, positions)
    for i in 1:N
        for j in i+1:N
            dists[i,j] == Inf && continue
            dist = _compute_interaction(interaction, geom, positions[:,i], positions[:,j])
            res[i,j] = dist
            res[j,i] = dist
        end
    end
    return res
end


### Implementation: ConstantInteraction
struct ConstantInteraction{T} <: Interaction
    val::T
end
_compute_interaction(c::ConstantInteraction, _, _, _) = c.val


### Implementation: PowerLaw
struct PowerLaw <: Interaction
    α::Float64
end
_compute_interaction(p::PowerLaw, geom, p1, p2) = _distance(geom, p1,p2)^(-p.α)


### Implementation: NNInteraction
struct NNInteraction{I} <: Interaction
    interaction::I
    k::Int
    NNInteraction(interaction; k=1) = new{typeof(interaction)}(interaction, k)
end

## Helpers
nearest_neighbour_from_interactions(interaction_mat, k=1) = nearest_neighbour_from_interactions!(copy(interaction_mat), k)

function nearest_neighbour_from_interactions!(interaction_mat, k=1)
    cutoffs = [sort!(unique(col); rev=true)[k] for col in eachcol(interaction_mat)] # k largest interactions
    interaction_mat[interaction_mat .< cutoffs'] .= 0 # not nearest nearest_neighbour = infinitely far away
    # could also consider converting to sparse but probably not worth it most of the times?
    # for sparse matrices simply remove the entry and not set it to Inf?
    return interaction_mat
end

## overwrite default
function interaction_matrix(nn::NNInteraction, geom, positions)
    J = interaction_matrix(nn.interaction, geom, positions)
    return nearest_neighbour_from_interactions!(J, nn.k)
end
