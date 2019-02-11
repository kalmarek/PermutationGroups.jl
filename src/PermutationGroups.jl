module PermutationGroups

using AbstractAlgebra
using Markdown

export AbstractOrbit, Orbit, Transversal, Schreier, StabilizerChain, PrmGroup
export firstmoved, fixes, fixedpoints
export base, getinv, representative, schreier_sims, sgs, sift

include("types.jl")
include("utils.jl")
include("orbit.jl")
include("words.jl")
include("schreier.jl")
include("stabchain.jl")
include("schreier-sims.jl")
include("groups.jl")

end # module
