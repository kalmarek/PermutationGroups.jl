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
Base.push!(esd::EigenSpaceDecomposition, a...) = push!(esd.eigenspaces, a...)
Base.append!(esd::EigenSpaceDecomposition, a) = append!(esd.eigenspaces)
LinearAlgebra.isdiag(esd::EigenSpaceDecomposition) = all(es -> isone(size(es,2)), esd)

function _restrict(M::MatrixElem{R}, basis::MatrixElem{R}) where R <: RingElem
    return basis'*M*basis
end

_dim(esd::EigenSpaceDecomposition) = isempty(esd.eigenspaces) ? 0 : sum(_dim, esd)
_dim(M::MatrixElem) = size(M,2)

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
