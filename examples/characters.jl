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
    eigen = Dict{elem_type(F), Generic.MatSpaceElem{GF}}()
    ct = 0
    for i in 0:order(F)-1
        if ct < size(M, 1)
            e = F(i)
            nullity, basis = nullspace(M - e*Id)
            if nullity > 0
                ct += nullity
                eigen[e] = basis
            end
        else
            break
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



struct EigenSpaceDecomposition{GF <: FinFieldElem}
    eigenspaces::Vector{Generic.MatSpaceElem{GF}}
end
Base.getindex(esd::EigenSpaceDecomposition, s) = getindex(esd.eigenspaces, s)
Base.iterate(esd::EigenSpaceDecomposition) = iterate(esd.eigenspaces)
Base.iterate(esd::EigenSpaceDecomposition, s) = iterate(esd.eigenspaces, s)
Base.length(esd::EigenSpaceDecomposition) = length(esd.eigenspaces)
_is_diagonal(esd::EigenSpaceDecomposition) = all((size.(esd, 2)) .== 1)

function eigen_space_decomposition(M::Generic.MatSpaceElem{GF}) where GF <: FinFieldElem
    return EigenSpaceDecomposition(collect(values(eigen(M))))
end

function _restrict(M::Generic.MatSpaceElem{GF}, basis) where GF <: FinFieldElem
    return basis*M*basis'
end

function eigen_space_decomposition(
                                  M::Generic.MatSpaceElem{GF}, 
                                  esd::EigenSpaceDecomposition{GF}) where GF <: FinFieldElem
    esdvector = typeof(esd.eigenspaces)()
    for basis in esd
        if size(basis, 2) == 1
            push!(esdvector, basis)
        else
            nesd = eigen_space_decomposition(_restrict(M, basis))
            for nbasis in nesd
                push!(esdvector, basis*nbasis)
            end
        end
        return EigenSpaceDecomposition(esdvector)
    end
end

function sd_basis(Ns::Vector{PermutationGroups.CCMatrix{T, C}}, F::GF) where {T, C, GF <: FinField}
    @assert !isempty(Ns)
    esd = eigen_space_decomposition(matrix(F, Ns[1]))
    ct = 1
    while !_is_diagonal(esd)
        ct += 1
        esd = eigen_space_decomposition(matrix(F, Ns[ct]), esd)
        if ct > length(Ns)
            @error "Matrices are not simultaneously diagonalizable"
        end
    end
    return esd
end

function _esd2matrix(esd::EigenSpaceDecomposition) 
    
end

basis = sd_basis(Ns, F)

#=
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
=#


