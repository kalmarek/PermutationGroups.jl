using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using AbstractAlgebra
using Revise
using PermutationGroups


# include("../src/characters.jl")
# RG = GroupRing(G)
# Build group ring for the permutation group over 4 elements.

N = 4
G = SymmetricGroup(N)
S = gens(G)
ccG = conjugacy_classes(G)

# Multiplication tables for conjugacy classes.
Ns = [PermutationGroups.CCMatrix(ccG, i) for i in 1:length(ccG)]

using LinearAlgebra

function LinearAlgebra.eigvals(M::Generic.MatSpaceElem{GF}) where GF<:FinFieldElem
    F = base_ring(M)
    Id = identity_matrix(M)
    eigvals = Dict{elem_type(F), Int}()
    for i in 0:order(F)-1
        e = F(i)
        nullity, basis = nullspace(M - e*Id)
        if nullity > 0
            eigvals[e] = nullity
        end
        # early break
    end
    return eigvals
end

function LinearAlgebra.eigen(M::Generic.MatSpaceElem{GF}) where GF<:FinFieldElem
    F = base_ring(M)
    Id = identity_matrix(M)
    eigen = Dict{elem_type(F), Any}()
    for i in 0:order(F)-1
        if length(eigen) < size(M, 1)
            e = F(i)
            nullity, basis = nullspace(M - e*Id)
            if nullity > 0
                eigen[e] = basis
            end
        end
    end
    return eigen 
end

F = GF(PermutationGroups.dixon_prime(ccG))
eigvals(matrix(F, Ns[1]))
eigvals(matrix(F, Ns[2]))
eigvals(matrix(F, Ns[3]))
eigvals(matrix(F, Ns[4]))
eigvals(matrix(F, Ns[5]))

eigs = let p = 23, Ms = M
    F = AbstractAlgebra.GF(p)
    AAMs = [matrix(F, m) for m in Ms] # if we implement our own nullspace we don't need this
    Id = identity_matrix(first(AAMs))

    potential_eigenvalues = Set(F(i) for i in 0:p-1)
    eigenvectors = Dict{elem_type(F), Vector{elem_type(F)}}()

    for m in AAMs[5:5]
        broken = false
        for e in potential_eigenvalues
            nullity, basis = nullspace(m - e*Id)
            if nullity > 0
                @show e, nullity, basis
                # find the "best" eigenvector among basis columns and add it to eigenvectors
                # basis.entries is a standard julia Array
                # entries = basis.entries[:,end]
                eigenvectors[e] = entries
                if length(eigenvectors) == size(m,1)
                    broken = true
                    break
                end
                delete!(potential_eigenvalues, e)
            end
        end
        broken && break
    end

    eigenvectors
end
