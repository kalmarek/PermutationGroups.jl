@static if VERSION < v"1.7"
    mutable struct Perm{T<:Integer} <: AbstractPermutation
        images::Vector{T}
        inv::Perm{T}
        cycles::CycleDecomposition{T}

        function Perm{T}(
            images::AbstractVector{<:Integer},
            check::Bool = true,
        ) where {T}
            check &&
                @assert isperm(images) "Provided images are not permutation!"
            deg = __degree(images)
            if deg == length(images)
                return new{T}(images)
            else
                # for future: use @time one(Perm{Int})
                # to check if julia can elide the creation of view
                return new{T}(@view images[Base.OneTo(deg)])
            end
        end
    end
else
    mutable struct Perm{T<:Integer} <: AbstractPermutation
        images::Vector{T}
        @atomic inv::Perm{T}
        @atomic cycles::CycleDecomposition{T}

        function Perm{T}(
            images::AbstractVector{<:Integer},
            check::Bool = true,
        ) where {T}
            check &&
                @assert isperm(images) "Provided images are not permutation!"
            deg = __degree(images)
            if deg == length(images)
                return new{T}(images)
            else
                # for future: use @time one(Perm{Int})
                # to check if julia can elide the creation of view
                return new{T}(@view images[Base.OneTo(deg)])
            end
        end
    end
end

# performance ?
@static if VERSION < v"1.7"
    function Base.copy(p::Perm)
        imgs = copy(p.images)
        q = typeof(p)(imgs, false)
        if isdefined(p, :inv)
            inv_imgs = copy(p.inv.images)
            q⁻¹ = typeof(p)(inv_imgs, false)
            q.inv = q⁻¹
            q⁻¹.inv = q
        end
        return q
    end

else
    function Base.copy(p::Perm)
        imgs = copy(p.images)
        q = typeof(p)(imgs, false)
        if isdefined(p, :inv, :sequentially_consistent)
            inv_imgs = copy(@atomic(p.inv).images)
            q⁻¹ = typeof(p)(inv_imgs, false)
            @atomic q.inv = q⁻¹
            @atomic q⁻¹.inv = q
        end
        return q
    end
end

# convienience constructor: inttype(::Type{<:AbstractPermutation}) defaults to UInt32
function Perm(images::AbstractVector{<:Integer}, check = true)
    return Perm{inttype(Perm)}(images, check)
end

# convienience conversion:
function Base.convert(::Type{Perm{T}}, p::Perm) where {T}
    inttype(p) == T && return p
    return Perm{T}(convert(Vector{T}, p.images), false)
end

# inttype must be T (UInt16 by default) since we store it in e.g. cycles
inttype(::Type{Perm{T}}) where {T} = T

# ## Interface of AbstractPermutation
degree(σ::Perm) = length(σ.images)

function Base.:^(n::Integer, σ::Perm)
    return 1 ≤ n ≤ degree(σ) ? oftype(n, @inbounds σ.images[n]) : n
end
@static if VERSION < v"1.7"
    function Base.inv(σ::Perm)
        if !isdefined(σ, :inv)
            σ⁻¹ = typeof(σ)(invperm(σ.images), false)
            σ.inv = σ⁻¹
            σ⁻¹.inv = σ
        end
        return σ.inv
    end
else
    function Base.inv(σ::Perm)
        if !isdefined(σ, :inv, :sequentially_consistent)
            σ⁻¹ = typeof(σ)(invperm(σ.images), false)
            # we don't want to end up with two copies of inverse σ floating around
            if !isdefined(σ, :inv, :sequentially_consistent)
                @atomic σ.inv = σ⁻¹
                @atomic σ⁻¹.inv = σ
            end
        end
        return σ.inv
    end
end

# optional
@static if VERSION < v"1.7"
    function cycles(σ::Perm)
        if !isdefined(σ, :cycles)
            cdec = CycleDecomposition(σ)
            σ.cycles = cdec
        end
        return σ.cycles
    end
else
    function cycles(σ::Perm)
        if !isdefined(σ, :cycles, :sequentially_consistent)
            cdec = CycleDecomposition(σ)
            # we can afford producing more than one cycle decomposition
            @atomic σ.cycles = cdec
        end
        return σ.cycles
    end
end

function parity(σ::Perm)
    isdefined(σ, :cycles) && return parity(cycles(σ))
    return __parity_generic(σ)
end
