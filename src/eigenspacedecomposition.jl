function LinearAlgebra.eigen(M::Generic.MatrixElem{GF}) where GF <: FinFieldElem
    @assert ==(size(M)...)
    F = base_ring(M)
    Id = identity_matrix(M)
    eigen = Dict{Int, typeof(M)}()
    cumdim = 0
    for i in 0:Int(order(F))-1
        cumdim >= _dim(M) && break
        nullity, basis = nullspace(M - i*Id) # looking for right eigenspaces
        if nullity > 0
            @debug "eigenvalue over F found: " e nullity, basis
            cumdim += nullity

            # for j in 1:ncols(basis)
            j = 1
            for i in 1:nrows(basis)
                if !iszero(basis[i, j]) && !isone(basis[i, j])
                    multiply_column!(basis, inv(basis[i,j]), j)
                    break
                end
            end
            # end

            eigen[i] = basis
        end
    end
    return eigen
end

_dim(M::MatrixElem) = size(M, 2) # we're right eigenspaces, hence column-based
_dims_to_ptrs!(dims) = (pushfirst!(dims, 1); cumsum!(dims, dims))

function eigen_decomposition!(M::MatrixElem{R}, eigspace_ptrs=Vector{Int}()) where R <: RingElement
    eigspace_ptrs = Vector{Int}()
    eigenvalues = eigen(M)
    sizehint!(eigspace_ptrs, length(eigenvalues)+1)
    push!(eigspace_ptrs, 1)
    for basis in values(eigenvalues)
        dim = _dim(basis)
        cd = eigspace_ptrs[end]
        ran = cd:cd+dim-1
        M[:, ran] = basis # we're column based
        push!(eigspace_ptrs, cd+dim)
    end
    @assert eigspace_ptrs[end] == _dim(M)+1
    return M, eigspace_ptrs
end

EigenSpaceDecomposition(R::Ring, nrows::Integer, ncols::Integer) =
    new{typeof(R)}(identity_matrix(R, dim, dim), [1, ncols+1])

EigenSpaceDecomposition(M::MatrixElem) =
    EigenSpaceDecomposition(eigen_decomposition!(deepcopy(M))...)

function Base.show(io::IO, ::MIME"text/plain", esd::EigenSpaceDecomposition)
    R = parent(esd.basis)
    println(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", parent(esd.basis))
    print(io, esd.basis)
end

function Base.show(io::IO, esd::EigenSpaceDecomposition{R}) where R
    print(io, tuple(diff(esd.eigspace_ptrs)...), "-splitting over ", parent(esd.basis))
end

Base.length(esd::EigenSpaceDecomposition) = length(esd.eigspace_ptrs)-1

function getindex(esd::EigenSpaceDecomposition, i::Integer)
    @boundscheck 1<= i <= length(esd)
    return esd.basis[:, esd.eigspace_ptrs[i]:esd.eigspace_ptrs[i+1]]
end

function Base.iterate(esd::EigenSpaceDecomposition, s=1)
    s > length(esd) && return nothing
    first_last = esd.eigspace_ptrs[s]:esd.eigspace_ptrs[s+1]-1
    return (esd.basis[:, first_last], s+1)
end

# right eigenspaces for left, (row-)eigenspaces), inv switches places
_change_basis(M::MatrixElem, basis::MatrixElem) = inv(basis)*M*basis

Base.eltype(::EigenSpaceDecomposition{R, M}) where {R, M} = M

LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) =
    all(isone, diff(esd.eigspace_ptrs))

function refine!(esd::EigenSpaceDecomposition{R}, M::MatrixElem{R}) where R
    m = _change_basis(M, esd.basis)
    @debug "matrices: original (M) and after base change (m)" M m
    m, eigspace_ptrs = eigen_decomposition!(m, esd.eigspace_ptrs)

    @assert size(m) == size(esd.basis)
    for i in 1:nrows(esd.basis)
        for j in 1:ncols(esd.basis)
            esd.basis[i,j] = m[i,j]
        end
    end

    resize!(esd.eigspace_ptrs, length(eigspace_ptrs))
    esd.eigspace_ptrs .= eigspace_ptrs
    return esd
end

function sd_basis(Ns, ring::Ring)
    @assert !isempty(Ns)
    esd = EigenSpaceDecomposition(matrix(ring, first(Ns)))
    for N in Iterators.rest(Ns, 2)
        esd = refine!(esd, matrix(ring, N))
        @debug N esd.eigspace_ptrs
        isdiag(esd) && return esd
    end
    return esd
end
