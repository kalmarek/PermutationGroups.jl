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

function EigenSpaceDecomposition(M::MatrixElem{R}) where R <: RingElement
    esd = EigenSpaceDecomposition(collect(values(eigen(M))))
    _dim(esd) == _dim(M) && return esd
    @warn "the subspace of dimension $(_dim(M)) does not fully split over $(R)" dims=(_dim.(esd))
    return EigenSpaceDecomposition([M])
end

function Base.show(io::IO, esd::EigenSpaceDecomposition)
    print(io, "Splitting of a vector space of dimension ", _dim(esd), "into")
    println(io, _dim.(esd), "-dimensional subspaces:" )
    for subspace in esd
        println(io, subspace)
    end
end

Base.getindex(esd::EigenSpaceDecomposition, s) = getindex(esd.eigenspaces, s)
Base.iterate(esd::EigenSpaceDecomposition) = iterate(esd.eigenspaces)
Base.iterate(esd::EigenSpaceDecomposition, s) = iterate(esd.eigenspaces, s)
Base.length(esd::EigenSpaceDecomposition) = length(esd.eigenspaces)
LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) = all(es -> isone(size(es,2)), esd)

function _restrict(M::MatrixElem{R}, basis::MatrixElem{R}) where R <: RingElem
    return basis'*M*basis
end

_dim(esd::EigenSpaceDecomposition) = isempty(esd.eigenspaces) ? 0 : sum(_dim, esd)
_dim(M::MatrixElem) = size(M,2)

function eigen_space_decomposition(
                                  M::Generic.MatSpaceElem{GF},
                                  esd::EigenSpaceDecomposition{GF}) where GF <: FinFieldElem
    eigenspaces = EigenSpaceDecomposition(Vector{typeof(M)}())
    sizehint!(eigenspaces, length(esd))
    dim = _dim(esd)
    @debug esd
    for subspace in esd
        if _dim(subspace) == 1
            push!(eigenspaces, subspace)
        else
            nesd = EigenSpaceDecomposition(_restrict(M, subspace))
            @debug _dim.(nesd)
            if _dim(nesd) == _dim(subspace) #length(nesd) > 0
                for nsubspace in nesd
                    push!(eigenspaces, subspace*nsubspace)
                end
            else
                @warn "the subspace does not fully split!: $(_dim.(nesd))"
                push!(eigenspaces, subspace)
            end
        end
    end
    nesd = EigenSpaceDecomposition(eigenspaces)
    @debug _dim.(eigenspaces)

    @assert length(nesd) >= length(esd)
    @assert _dim(nesd) == dim
    return nesd
end

function sd_basis(Ns::Vector{CCMatrix{T, C}}, ring::Ring) where {T, C}
    @assert !isempty(Ns)
    esd = EigenSpaceDecomposition([matrix(F, Ns[1])])
    ct = 1
    while !isdiag(esd)
        ct += 1
        esd = eigen_space_decomposition(matrix(F, Ns[ct]), esd)
        if ct > length(Ns)
            @error "Matrices are not simultaneously diagonalizable"
        end
        @info [isone(size(es,2)) for es in esd]
    end
    return esd
end
