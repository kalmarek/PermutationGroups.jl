module PermutationGroups

using AbstractAlgebra
using Markdown

import AbstractAlgebra: AbstractPermutationGroup, AbstractPerm, Group, GroupElem, Generic.Perm,
    mul!, @perm_str, degree, order, permtype, gens

export AbstractPermutationGroup, AbstractPerm, Group, GroupElem
export @perm_str, Perm, degree, order, permtype, gens

export AbstractOrbit, Orbit, Transversal, Schreier, StabilizerChain, PermGroup
export firstmoved, fixes, fixedpoints, lastmoved
export base, getinv, representative, schreier_sims, sgs, sift

include("types.jl")
include("utils.jl")
include("orbit.jl")
include("words.jl")
include("schreier.jl")
include("groups.jl")
include("stabchain.jl")
include("schreier-sims.jl")

end # module
