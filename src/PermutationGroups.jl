module PermutationGroups

using GroupsCore
using Markdown

import GroupsCore.AbstractAlgebra:
    AbstractPermutationGroup,
    AbstractPerm,
    Generic.Perm,
    @perm_str,
    permtype

export @perm_str, Perm, degree, order, permtype, gens

export AbstractOrbit, Orbit, Transversal, Schreier, StabilizerChain, PermGroup
export firstmoved, fixes, fixedpoints, lastmoved, nfixedpoints
export base, representative, schreier_sims, sgs, sift

include("types.jl")
include("orbit.jl")
include("words.jl")
include("group_interface.jl")
include("groups.jl")
include("schreier.jl")
include("stabchain.jl")
include("schreier-sims.jl")

include("utils.jl")

end # module
