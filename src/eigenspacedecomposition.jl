function LinearAlgebra.eigen(M::Generic.MatSpaceElem{GF}) where GF <: FinFieldElem
    @assert ==(size(M)...)
    F = base_ring(M)
    Id = identity_matrix(M)
    eigen = Dict{elem_type(F), typeof(M)}()
    count = 0
    for i in 0:order(F)-1
        count >= size(M, 1) && break

        e = F(i)
        nullity, basis = nullspace(M - e*Id)
        @debug e nullity, basis
        if nullity > 0
            count += nullity
            !iszero(basis[1,1]) && multiply_column!(basis, inv(basis[1,1]), 1)
            eigen[e] = basis
        end
    end
    return eigen
end

_dim(M::MatrixElem) = size(M, 2) # we're column based
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
        M[:, ran] = basis
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
    R = parent(first(esd))
    println(io, diff(esd.eigspace_ptrs),
        " - spliting of $R module of dim $(_dim(esd.basis)).")
    println(io, esd.basis)
end

function Base.show(io::IO, esd::EigenSpaceDecomposition{R}) where R
    print(io, diff(esd.eigspace_ptrs), "-eigenspace splitting over ", R)
end

Base.length(esd::EigenSpaceDecomposition) = length(esd.eigspace_ptrs)-1
_change_basis(M::MatrixElem, basis::MatrixElem) = inv(basis)*M*basis

LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) =
    all(isone, diff(esd.eigspace_ptrs))

function refine!(esd::EigenSpaceDecomposition{R}, M::MatrixElem{R}) where R
    m = _change_basis(M, esd.basis)
    @debug "matrices: original (M) and after base change (m)" M m
    m, eigspace_ptrs = eigen_decomposition!(m, esd.eigspace_ptrs)
    esd.basis.entries .= m.entries
    resize!(esd.eigspace_ptrs, length(eigspace_ptrs))
    esd.eigspace_ptrs .= eigspace_ptrs
    return esd
end

function sd_basis(Ns::Vector{CCMatrix{T, C}}, ring::Ring) where {T, C}
    @assert !isempty(Ns)
    esd = EigenSpaceDecomposition([matrix(ring, first(Ns))])

    for N in Iterators.rest(Ns, 2)
        esd = refine(esd, matrix(ring, N))
        @debug N [isone(size(es,2)) for es in esd]
        isdiag(esd) && return esd
    end
    return esd
end
