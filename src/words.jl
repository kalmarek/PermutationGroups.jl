################################################################
#
# representation of words as vector of pointers to generators
#
################################################################

@inline Base.push!(pw::Word, i) = push!(pw.gens_ptrs, i)
@inline Base.getindex(pw::Word, n) = pw.gens_ptrs[n]
@inline Base.setindex!(pw::Word, val, n) = pw.gens_ptrs[n] = val

@inline Base.iterate(pw::Word) = iterate(pw.gens_ptrs)
@inline Base.iterate(pw::Word, s) = iterate(pw.gens_ptrs, s)
@inline Base.length(pw::Word) = length(pw.gens_ptrs)
@inline Base.eltype(::Word{I}) where I = I

@inline function Base.:(==)(pw1::Word, pw2::Word)
    return pw1.gens_ptrs == pw2.gens_ptrs
end

Base.hash(pw::Word, h::UInt) = hash(Word, hash(pw.gens_ptrs, h))

@doc doc"""
    Word(orb, gens_inv, pt, [^])
> return a `Word` `w` consising of pointers to inverses of generators `gens_inv`.
> The group element `r` such that `first(orb)^r = pt` can be
> reconstructed via `prod(gens_inv[i] for i in reverse(w))`
> e.g. if
    g = gens[1]*gens[3]*gens[2]
    δ = first(orb)^g
> then the result _may be_ `[2,3,1]`.
"""
function Word(gens_inv::Vector{<:GroupElem}, orb::AbstractOrbit{T, I}, pt, op=^) where {T, I<:Integer}
    γ = pt
    res_ptr = I[]
    while γ ≠ first(orb)
        idx = orb[γ]
        push!(res_ptr, idx)
        γ = op(γ, gens_inv[idx]) #assume that this is already the inverse!
    end
    return Word(res_ptr)
end

@doc doc"""
    (w::Word)(gens, [init])
> Computes the evaluation of group word `w` as a group element in generators
> `gens`. Optional `init` element is by default the identity.
"""
function (pw::Word)(gens::Vector{<:GroupElem}, init=Perm(degree(first(gens))))
    res = init
    @inbounds for i in pw
        res = fastmul!(res, res, gens[i])
    end
    return res
end

@doc doc"""
    representative(orb, gens, pt, [^])
> Return a representative of the coset of `Stab(orb)` which takes `first(orb)`
> to `pt` i.e. an element `g` of `⟨gens⟩` such that `first(orb)^g = pt`.
"""
function representative(gens::Vector{<:GroupElem}, orb::AbstractOrbit{I, <:Integer}, pt::I, op=^) where I<:Integer
    gens_inv = inv.(gens)
    perm_word = Word(gens_inv, orb, pt, op)
    g = perm_word(gens_inv)
    return inv(g)
end 
