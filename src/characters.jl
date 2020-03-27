using LinearAlgebra
using AbstractAlgebra

#=
using GroupRings
function Base.getindex(M::CCMatrix{T, C}, s::Integer, t::Integer) where {T, C <: GroupRingElem}
    if isone(-M.m[s,t])
        r = M.r
        res = M.cc[r]*M.cc[s]
        for (t, cc) in enumerate(M.cc)
            g = first(supp(cc)) # it doesn't matter which we take
            M.m[s,t] = res[g] # c_RST = r*s = t
        end
    end
    return M.m[s,t]
end
=#

using DynamicPolynomials
using MultivariatePolynomials

function characteristic_polynomial(M::AbstractArray, T::MultivariatePolynomials.AbstractVariable)
    n,m = size(M)
    @assert n==m
    L = convert.(polynomialtype(T), -M)
    for i in 1:n
        L[i,i]+=T
    end
    return det(L)
end

function Base.mod(p::MultivariatePolynomials.AbstractPolynomialLike, m::Int)
    return dot(mod.(coefficients(p), m), MultivariatePolynomials.monomials(p))
end

function eval_mod(p::MultivariatePolynomials.AbstractPolynomialLike, m::Int)
    return [mod(p(i), m) for i = 0:m-1]
end

function factor_mod(p::MultivariatePolynomials.AbstractPolynomialLike, m::Int)
    values = eval_mod(p, m)
    println(values)
    if all(values.!=0)
        return p
    else
        T = first(variables(p))
        id = findall(x->x==0, values)
        factors = [T-mod(i-1,m) for i in id]
        println(factors)
        r = p
        for f in factors
            r,_ = divrem(r,f)
            println(r)
        end
    end
    println("here")
    return [convert.(typeof(p),factors)..., factor_mod(convert(typeof(p),r),m)]
end

"""

eigenvectors = Dict{Polynomial, Vector}()[] # pseudocode
for cc in ccG
  M = CCMat(cc, i) # M is class multiplication matrix 
  F ← minimal polynomial of M
  for f in factor(F) # factors correspond to distinct eigenvalues/eigenvectors
    haskey(eigvecs, f) && continue
    K = nullspace(f(M)) # kernel corresponds to eigenspace
    v = ... # use QR/lll to find "best" integer vector in K
    eigenvectors[f] = v
    length(eigvecs) == size(M, 1) && break
  end
  length(eigvecs) == size(M, 1) && break
  # double break could be avoided by using a flag or so
end

then we need to "normalize" eigenvectors and we're done;

* NOTE: if p is small it's cheaper to keep the list of eigenvalues instead of
factors (and instead computing F/its factorisation) ask for the nullspace of
M-eI where e is a new eigenvalue. Maybe this is even better in general *

To lift them to C (or R if they're all real) we need to do a little bit more,
but not much ;)
"""

function LinearAlgebra.dot(a::AbstractVector{T}, b::AbstractVector{T}) where T <: AbstractAlgebra.GFElem
    @assert !isempty(a)
    @assert length(a) == length(b)
    return sum(a[i]*b[i] for i = 1:length(a))
end

"""
    _projection(a, u) 

Returns projection of a into the subspace spanned by u.
"""
function _projection(a::AbstractVector, u::AbstractVector)
    return dot(a, u)*inv(dot(u, u)).*u
end

"""
    _normalize(a::AbstractVector)

Returns a normalized to a unit vector.
"""
function _normalize(a::AbstractVector)
    l = dot(a, a)
    @assert l!=0
    return a.*inv(l)
end

"""
    gram_schmidt(A::AbstractMatrix)

Return matrices Q and R such that A = QR, where Q is orthogonal and R is an upper triangular Matrix, using the Gram-Schmidt process.
"""
function gram_schmidt(A::AbstractMatrix)
    n, m = size(A)
    @assert n == m 
    Q = zeros(parent(A[1,1]), n, n)
    R = zeros(parent(A[1,1]), n, n)

    Q[:, 1] = _normalize(A[:, 1])
    R[1, 1] = dot(Q[:, 1], A[:, 1])

    @info Q[:, 1]

    for k = 2:n
   #     Q[:, k] = _normalize(A[:, k] - sum(_projection(A[:, k], Q[:, j]) for j = 1:k-1))
        Q[:, k] = A[:, k] - sum(_projection(A[:, k], Q[:, j]) for j = 1:k-1)
        @info Q[:, k]
        @info dot(Q[:, k], Q[:, k])
        for j = 1:k
            R[j, k] = dot(Q[:, j], A[:, k])
        end
    end
    return Q, R
end

