mutable struct Perm{T<:Integer} <: AbstractPermutation
    images::Vector{T}
    @atomic inv::Perm{T}
    @atomic cycles::CycleDecomposition{T}

    function Perm{T}(
        images::AbstractVector{<:Integer},
        check::Bool = true,
    ) where {T}
        if check
            isperm(images) || error("Provided images are not permutation!")
        end

        deg = __degree(images)
        σ = new{T}(@view images[firstindex(images):deg])

        return σ
    end
end

# convienience constructor: defaults to UInt32
Perm(images::AbstractVector{<:Integer}) = Perm{inttype(Perm)}(images)

# inttype must be T (UInt32 by default) since we store it in e.g. cycles
inttype(::Type{Perm{T}}) where {T} = T

# ## Interface of AbstractPermutation
degree(σ::Perm) = length(σ.images)

function Base.:^(n::Integer, σ::Perm)
    return 1 ≤ n ≤ degree(σ) ? oftype(n, @inbounds σ.images[n]) : n
end

function Base.inv(σ::P) where {P<:Perm}
    # if @atomic !isdefined(σ, :inv)
    if !isdefined(σ, :inv)
        σ_inv = P(invperm(σ.images), false)
        # we don't want to end up with two copies of inverse σ floating around
        # if @atomic !isdefined(σ, :inv)
        if !isdefined(σ, :inv)
            @atomic σ.inv = σ_inv
            @atomic σ.inv.inv = σ
        end
    end
    return σ.inv
end

# optional

function cycles(σ::Perm)
    # if @atomic !isdefined(σ, :cycles)
    if !isdefined(σ, :cycles)
        cdec = CycleDecomposition(σ)
        # we can afford producing more than one cycle decomposition
        @atomic σ.cycles = cdec
    end
    return σ.cycles
end
