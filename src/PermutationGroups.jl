module PermutationGroups

using GroupsCore
using Random
using Markdown

include("Perms/Perms.jl")
import .Perms:
    Perm, AbstractPermutation, degree, firstmoved, permtype, @perm_str

export @perm_str, Perm, permtype

export AbstractOrbit,
    Orbit, PermGroup, Permutation, Schreier, StabilizerChain, Transversal
export firstmoved, fixes, fixedpoints, lastmoved, nfixedpoints
export base, degree, representative, schreier_sims, sgs, sift

abstract type AbstractPermutationGroup <: Group end

include("types.jl")
include("orbit.jl")
include("words.jl")
include("io.jl")
include("group_interface.jl")
include("groups.jl")
include("schreier.jl")
include("stabchain.jl")
include("schreier-sims.jl")

include("actions.jl")
end # module
