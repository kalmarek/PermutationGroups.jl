module Perms

using GroupsCore
export AbstractPermutation, Perm, cycles, degree, permtype, @perm_str

include("cycle_decomposition.jl")

# abstract definitions
include("abstract_perm.jl")
include("arithmetic.jl")
include("misc.jl")

include("utils.jl")

# concrete implementations
include("perm_images.jl")
include("macro_perm.jl")
end
