"""
    AbstractPermutation
Abstract type representing permutations of set `1:n`.

Subtypes `Perm <: AbstractPermutation` must implement the following functions:
* `Base.:^(i::Integer, σ::Perm)` - the image of `i` under `σ`,
* `degree(σ::Perm)` - the minimal `n` such that `k^σ == k` for all `k > n`,
* `Perm(images::AbstractVector{<:Integer}[, check::Bool=true])` - construct a
`Perm` from a vector of images. Optionally the second argument `check` may be
set to `false` when the caller knows that `images` constitute a honest
permutation.

!!! Note
    There is no formal requirement that the `Perm(images)` constructor actually
    returns a `Perm`. Any `AbstractPermutation` object would do. This may be
    useful if constructing permutation from images is not technically feasible.
"""
abstract type AbstractPermutation <: GroupElement end

"""
    degree(σ::AbstractPermutation)
Return a minimal number `n` such that `σ(k) == k` for all `k > n`.

Such number `n` can be understood as a _degree_ of a permutation, since we can
regard `σ` as an element of `Sym(n)` (and not of `Sym(n-1)`).

By convention `degree` of the trivial permutation must return `1`.
"""
function degree(σ::AbstractPermutation)
    throw("not implemented: `PermutationGroups.degree(::$(typeof(σ)))`")
end

"""
    ^(i::Integer, σ::AbstractPermutation)
Return the image of `i` under the (permutation) action of `σ`.

We consider `σ` as a finite support permutation of `ℕ`, so by convention `k^σ = k`
for all `k > degree(σ)`.
"""
function Base.:^(::Integer, σ::AbstractPermutation)
    throw("not implemented: Base.:^(::Integer, ::$(typeof(σ)))")
end

Base.one(::Type{P}) where {P<:AbstractPermutation} = P(inttype(P)[1], false)
Base.one(σ::AbstractPermutation) = one(typeof(σ))
Base.isone(σ::AbstractPermutation) = degree(σ) == 1

inttype(::Type{P}) where {P<:AbstractPermutation} = UInt32
inttype(σ::AbstractPermutation) = inttype(typeof(σ))

function Base.:(==)(σ::AbstractPermutation, τ::AbstractPermutation)
    degree(σ) ≠ degree(τ) && return false
    for i in Base.OneTo(degree(σ))
        if i^σ != i^τ
            return false
        end
    end
    return true
end

function Base.hash(σ::AbstractPermutation, h::UInt)
    h = hash(AbstractPermutation, h)
    for i in Base.OneTo(degree(σ))
        h = hash(i^σ, h)
    end
    return h
end

Base.broadcastable(p::AbstractPermutation) = Ref(p)

cycles(σ::AbstractPermutation) = CycleDecomposition(σ)

function CycleDecomposition(σ::AbstractPermutation)
    T = inttype(σ)
    deg = degree(σ)

    # allocate vectors of the expected size
    visited = falses(deg)
    cycles = Vector{T}(undef, deg)
    # expected number of cycles - (overestimation of) the harmonic
    cyclesptr = Vector{T}(undef, 5 + ceil(Int, Base.log(deg + 1)))

    # shrink them accordingly
    resize!(cycles, 0)
    resize!(cyclesptr, 1)
    cyclesptr[begin] = 1

    @inbounds for idx in Base.OneTo(deg)
        visited[idx] && continue
        first_pt = idx

        orbit_len = 1
        push!(cycles, first_pt)
        visited[first_pt] = true
        next_pt = first_pt^σ
        while next_pt ≠ first_pt
            orbit_len += 1
            push!(cycles, next_pt)
            visited[next_pt] = true
            next_pt = next_pt^σ
        end
        push!(cyclesptr, last(cyclesptr) + orbit_len)
        # cycles[cyclesptr[i]:cyclesptr[i+1]-1] contains i-th cycle

        # Δ = orbit_plain(T(i), σ, ^)
        # visited[Δ] .= true
        # append!(cycles, Δ)
        # push!(cyclesptr, cyclesptr[end] + length(Δ))

    end
    return CycleDecomposition{T}(cycles, cyclesptr)
end

# function orbit_plain(x, s::GroupElement, action = ^)
#     Δ = [x]
#     γ = action(x, s)
#     while γ != x
#         push!(Δ, γ)
#         γ = action(γ, s)
#     end
#     return Δ
# end

# IO

Base.show(io::IO, g::AbstractPermutation) = _print_perm(io, g)

function _print_perm(
    io::IO,
    p::AbstractPermutation,
    width::Integer = last(displaysize(io)),
)
    if isone(p)
        return print(io, "()")
    else
        @assert width > 3
        cum_length = 0
        for c in cycles(p)
            length(c) == 1 && continue
            cyc = join(c, ",")

            if width - cum_length >= length(cyc) + 2
                print(io, "(", cyc, ")")
                cum_length += length(cyc) + 2
            else
                available = width - cum_length - 3
                print(io, "(", SubString(cyc, 1, available), " …")
                break
            end
        end
    end
end
