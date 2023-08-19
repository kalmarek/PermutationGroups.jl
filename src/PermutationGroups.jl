module PermutationGroups

using GroupsCore
import GroupsCore: istrivial
using Random

abstract type AbstractPermutationGroup <: Group end
Base.IteratorSize(::Type{<:AbstractPermutationGroup}) = Base.HasLength()

include("Perms/Perms.jl")
import .Perms:
    Perm, AbstractPermutation, cycles, degree, firstmoved, permtype, @perm_str

export @perm_str, Perm, permtype

export AbstractOrbit,
    Orbit,
    AbstractTransversal,
    Transversal,
    SchreierTransversal,
    PermGroup,
    Permutation,
    StabilizerChain
export firstmoved, fixes, fixedpoints, lastmoved, nfixedpoints
export base, degree, representative, schreier_sims, sgs, sift

include("orbit.jl")
include("stabchain.jl")
include("schreier_sims.jl")
include("perm_group.jl")
include("group_interface.jl")

include("precompile.jl")

end # module
