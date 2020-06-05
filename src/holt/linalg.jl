using AbstractAlgebra

export SubSpace
export depth_vector, vector_in_subspace, invert_matrix, nullspace
"""
    SubSpace{T}

Type to store a subspace. The basis Vectors are normalized such that 
the first nonzero entry of each vector is 1 and the index of the first 
non-zero entry is different for each vector in the basis. These positions 
are stored in the field positions. 
"""
mutable struct SubSpace{T}
    basis::Vector{Vector{T}}
    positions::Vector{Int}
end

_basis(W::SubSpace) = W.basis
_positions(W::SubSpace) = W.positions
SubSpace{T}() where T = SubSpace(Vector{T}[], Int[])


"""
    depth_vector(v::Vector{T})

Locates the leading nonzero coefficient α of a vector and returns 
its position and α itself, or 0, 0 if the vector is zero.
"""
function depth_vector(v::Vector{T}) where T
    @assert !isempty(v)
    α = zero(v[1])
    j = 0
    for i = 1:length(v)
        if !iszero(v[i])
            α = v[i]
            j = i
            break
        end
    end
    return j, α
end
"""
    vector_in_subspace(W::SubSpace{T}, v::Vector{T}; add = false)

Tests whether v is in W and returns true or false accordingly. If yes the coefficients
of v expressed in the basis of W are returned as a second argument.
If add = true, W is replaced by span(W,v) and a second output provides the coefficients 
of v expressed in the updated basis of W.
"""
function vector_in_subspace(W::SubSpace{T}, v::Vector{T}; add = false) where T
    β = _basis(W)
    γ = _positions(W)
    c = T[]
    for (i, p) in enumerate(γ)
        push!(c, v[p])
        iszero(v[p]) ? nothing : v -= v[p].*β[i]
    end
    d, a = depth_vector(v)
    if iszero(d)
        return true, c
    elseif add
        push!(β, inv(a).*v)
        push!(γ, d)
        push!(c, a)
        return false, c
    else
        return false
    end
end

function SubSpace(V::Vector{Vector{T}}) where {T}
    W = SubSpace{T}()
    for v in V
        vector_in_subspace(W, v; add = true)
    end
    return W
end

"""
    invert_matrix(A::Array{T, 2}) 

Returns the inverse of a square matrix A.
"""
function invert_matrix(AA::Array{T, 2}) where T
    if T <: AbstractFloat
        # This algorithm is numerically not stable, use default method instead.
        return AA^(-1)
    end
    @assert !isempty(AA)
    (d, l) = size(AA)
    @assert d==l "Matrix A is not square."
    I = one(AA[1,1]).*zeros(Int, d, d)
    for c = 1:d
        I[c, c] += 1
    end
    A = deepcopy(AA)
    for c = 1:d
        for s in [c, d]
            if iszero(A[s, c])
                throw(ErrorException("Matrix is not invertible."))
            end
            I[s, :] = inv(A[s, c]).*I[s, :]
            A[s, :] = inv(A[s, c]).*A[s, :]
            for t in 1:d
                if (t<c || t>s)
                    I[t, :] -= A[t, s].*I[s, :]
                    A[t, :] -= A[t, s].*A[s, :]
                end
            end
            if !(c==s)
                I[c, :], I[s, :] = I[s, :], I[c, :]
                A[c, :], A[s, :] = A[s, :], A[c, :]
            end
            break
        end
    end
    return I
end

"""
    nullspace(A)

Returns a SubSpace representing the Nullspace of A.
"""
function AbstractAlgebra.nullspace(A::Array{T, 2}) where T
    W = SubSpace{T}()
    B = transpose(deepcopy(A))
    (c, d) = size(B)
    r = 1
    l = Int[]
    for t = 1:d
        s = findfirst(x -> !iszero(x), B[r:c, t])
        if !(s isa Nothing)
            s += r-1
            B[s,:] = inv(B[s,t]).*B[s, :]
            for k = 1:c
                if !(k==s)
                    B[k,:] -= B[k,t].*B[s, :]
                end
            end
            if r!=s
                B[s,:], B[r,:] = B[r,:], B[s,:] 
            end
            push!(l, t) 
            r+=1 
        end
    end
    for t in setdiff(1:d, l)
        v = [zero(A[1,1]) for i = 1:d] 
        v[t] = one(A[1,1])
        for s = 1:r-1
            v[l[s]] = -B[s, t]
        end
        vector_in_subspace(W, v; add = true)
    end
    return W
end

