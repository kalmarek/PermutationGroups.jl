module PermutationGroups

using GroupsCore
import GroupsCore: istrivial
using Random
import AbstractPermutations as AP
import AbstractPermutations: AbstractPermutation

export Perm, SPerm, @perm_str, @sperm8_str

export AbstractOrbit,
    Orbit,
    Transversal,
    SchreierTransversal,
    PermGroup,
    Permutation,
    StabilizerChain
export base, representative, schreier_sims, sgs, sift

include("Perms/utils.jl")

# concrete implementations
include("Perms/perm_images.jl")
include("Perms/perm_static.jl")

abstract type AbstractPermutationGroup <: Group end
Base.IteratorSize(::Type{<:AbstractPermutationGroup}) = Base.HasLength()
Perms.degree(G::AbstractPermutationGroup) = maximum(degree, gens(G))

include("orbit.jl")
include("stabchain.jl")
include("schreier_sims.jl")
include("perm_group.jl")
include("group_interface.jl")
include("io.jl")

include("precompile.jl")

end # module
