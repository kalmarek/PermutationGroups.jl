module PermutationGroups

using GroupsCore
using Random
using Markdown

import GroupsCore.AbstractAlgebra:
    AbstractPermutationGroup,
    AbstractPerm,
    Generic.Perm,
    @perm_str,
    permtype

export @perm_str, Perm, permtype

export AbstractOrbit, Orbit, PermGroup, Permutation, Schreier, StabilizerChain, Transversal
export firstmoved, fixes, fixedpoints, lastmoved, nfixedpoints
export base, degree, representative, schreier_sims, sgs, sift

include("types.jl")
include("orbit.jl")
include("words.jl")
include("group_interface.jl")
include("groups.jl")
include("schreier.jl")
include("stabchain.jl")
include("schreier-sims.jl")

include("actions.jl")

end # module
