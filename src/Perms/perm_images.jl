@static if VERSION < v"1.7"
    mutable struct Perm{T<:Integer} <: AP.AbstractPermutation
        images::Vector{T}
        inv::Perm{T}
        cycles::CycleDecomposition{T}

        function Perm{T}(
            images::AbstractVector{<:Integer},
            check::Bool = true,
        ) where {T}
            if check && !isperm(images)
                throw(ArgumentError("Provided images are not permutation!"))
            end
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
    mutable struct Perm{T<:Integer} <: AP.AbstractPermutation
        images::Vector{T}
        @atomic inv::Perm{T}
        @atomic cycles::AP.CycleDecomposition{T}

        function Perm{T}(
            images::AbstractVector{<:Integer},
            check::Bool = true,
        ) where {T}
            if check && !isperm(images)
                throw(ArgumentError("Provided images are not permutation!"))
            end
            li = lastindex(images)
            deg = @inbounds images[li] ≠ li ? li : __degree(images)
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

# convienience constructor: inttype(::Type{<:AbstractPermutation}) defaults to UInt32
function Perm(images::AbstractVector{<:Integer}, check = true)
    return Perm{AP.inttype(Perm)}(images, check)
end

# ## Interface of AbstractPermutation
function Base.:^(n::Integer, σ::Perm)
    return 1 ≤ n ≤ AP.degree(σ) ? oftype(n, @inbounds σ.images[n]) : n
end
AP.degree(σ::Perm) = length(σ.images)
# inttype must be T (UInt16 by default) since we store it in e.g. cycles
AP.inttype(::Type{Perm{T}}) where {T} = T
AP.inttype(::Type{Perm}) = UInt16

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

    function Base.inv(σ::Perm)
        if !isdefined(σ, :inv)
            σ⁻¹ = typeof(σ)(invperm(σ.images), false)
            σ.inv = σ⁻¹
            σ⁻¹.inv = σ
        end
        return σ.inv
    end

    function AP.cycles(σ::Perm)
        if !isdefined(σ, :cycles)
            cdec = AP.CycleDecomposition(σ)
            σ.cycles = cdec
        end
        return σ.cycles
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

    function Base.inv(σ::Perm)
        if !isdefined(σ, :inv, :sequentially_consistent)
            if isone(σ)
                @atomic σ.inv = σ
            else
                σ⁻¹ = typeof(σ)(invperm(σ.images), false)
                # we don't want to end up with two copies of inverse σ floating around
                if !isdefined(σ, :inv, :sequentially_consistent)
                    @atomic σ.inv = σ⁻¹
                    @atomic σ⁻¹.inv = σ
                end
            end
        end
        return σ.inv
    end

    function AP.cycles(σ::Perm)
        if !isdefined(σ, :cycles, :sequentially_consistent)
            cdec = AP.CycleDecomposition(σ)
            # we can afford producing more than one cycle decomposition
            @atomic σ.cycles = cdec
        end
        return σ.cycles
    end
end

function Base.isodd(σ::Perm)
    isdefined(σ, :cycles) && return isodd(AP.cycles(σ))
    return AP.__isodd(σ)
end

macro perm_str(str)
    return :(Base.parse(Perm, $str))
end
