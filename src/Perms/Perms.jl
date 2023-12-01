module Perms

using GroupsCore
export AbstractPermutation, Perm, cycles, degree, permtype, @perm_str

include("utils.jl")

# concrete implementations
include("perm_images.jl")
include("macro_perm.jl")
end
