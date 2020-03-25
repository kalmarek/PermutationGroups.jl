using LinearAlgebra
using AbstractAlgebra
using GroupRings

struct CCMatrix{T, C} <: AbstractMatrix{T} # M_r
    cc::Vector{C} # vector of conjugacy classes to fix the order
    r::Int # the index of conjugacy class
    m::Matrix{T} # cache of class coefficients

    function CCMatrix(cc::A, r::Int, T::Type=Int) where {C, A<:AbstractVector{C}}
        M = -ones(T, length(cc), length(cc))
        new{T, C}(cc, r, M)
    end
end

Base.size(M::CCMatrix) = size(M.m)
Base.IndexStyle(::Type{<:CCMatrix}) = IndexCartesian()

function Base.getindex(M::CCMatrix{T, C}, s::Integer, t::Integer) where {T, C <: GroupRingElem}
    if isone(-M.m[s,t])
        r = M.r
        g = first(supp(M.cc[t])) # it doesn't matter which we take
        M.m[s,t] = (M.cc[r]*M.cc[s])[g] # c_RST = r*s = t
    # we compute much more above: no we obtain the whole row (? or column?)
    end
    return M.m[s,t]
end

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

function is_prime(p::Int)
    if p==1
        return false
    else
        prime = p==2 || p==3
        if !prime 
            prime = mod(p,6)==1||mod(p,6)==5
            if prime 
                for i = 5:floor(Int, sqrt(p),)
                    prime = mod(p, i)>0
                    if !prime
                        break
                    end
                end
            end
        end
        return prime
    end
end

function next_prime(z::Int)
    z += 1
    while !is_prime(z)
        z += 1
    end
    return z
end

function find_prime(m::Int) 
    @assert m < 1000000
    p = 2
    while mod(p, m)!=1
        p = next_prime(p)
    end
    return p
end

function find_int(m::Int, p::Int)
    z = 2
    while !(mod(z^m, p)==1) || any([mod(z^m, q)==1 for q in 1:p-1])
         z+=1
    end
    return z
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
