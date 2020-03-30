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

function _restrict(M::MatrixElem{R}, basis::MatrixElem{R}) where R <: RingElem
    return basis'*M*basis
end

LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) =
    all(isone, diff(esd.eigspace_ptrs))

function refine(esd::EigenSpaceDecomposition{R}, M::MatElem{R}) where R <: RingElem
    @debug esd

    new_esd = EigenSpaceDecomposition{R}()
    sizehint!(new_esd.eigenspaces, length(esd))

    for subspace in esd
        if isone(_dim(subspace))
            push!(new_esd, subspace)
        else
            sub_esd = EigenSpaceDecomposition(_restrict(M, subspace))
            @debug _dim.(sub_esd)
            for nsubspace in sub_esd
                push!(new_esd, subspace*nsubspace)
            end
        end
    end

    @debug _dim.(new_esd)

    @assert length(new_esd) >= length(esd)
    @assert _dim(new_esd) == _dim(esd)
    return new_esd
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
