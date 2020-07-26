module PermutationGroups

using LinearAlgebra
using AbstractAlgebra
using Primes
using Markdown

import AbstractAlgebra: AbstractPermutationGroup, AbstractPerm, Group, mul!

export AbstractOrbit, Orbit, Transversal, Schreier, StabilizerChain, PermGroup
export firstmoved, fixes, fixedpoints, lastmoved
export base, conjugacy_classes, getinv, representative, schreier_sims, sgs, sift

include("types.jl")
include("utils.jl")
include("orbit.jl")
include("words.jl")
include("schreier.jl")
include("groups.jl")
include("stabchain.jl")
include("schreier-sims.jl")

end # module
