mutable struct Perm{T<:Integer} <: AP.AbstractPermutation
    images::Vector{T}
    @atomic inv::Perm{T}
    @atomic cycles::AP.CycleDecomposition{T}

    function Perm{T}(images::Vector{T}; check::Bool = true) where {T}
        if check && !isperm(images)
            throw(ArgumentError("Provided images are not permutation!"))
        end
        li = lastindex(images)
        deg =
            iszero(li) ? li : @inbounds images[li] ≠ li ? li : __degree(images)
        resize!(images, deg)
        return new{T}(images)
    end
end

function Perm{T}(
    images::AbstractVector{<:Integer};
    check::Bool = true,
) where {T}
    deg = __degree(images)
    return Perm{T}(
        convert(Vector{T}, @view images[Base.OneTo(deg)]);
        check = check,
    )
end

# convienience constructor: inttype(::Type{<:AbstractPermutation}) defaults to UInt32
function Perm(images::AbstractVector{<:Integer}; check = true)
    return Perm{AP.inttype(Perm)}(images; check = check)
end

# ## Interface of AbstractPermutation
function Base.:^(n::Integer, σ::Perm)
    return 1 ≤ n ≤ AP.degree(σ) ? oftype(n, @inbounds σ.images[n]) : n
end
AP.degree(σ::Perm) = length(σ.images)
# inttype must be T (UInt16 by default) since we store it in e.g. cycles
AP.inttype(::Type{Perm{T}}) where {T} = T
AP.inttype(::Type{Perm}) = UInt16
AP.__unsafe_image(n::Integer, σ::Perm) = oftype(n, @inbounds σ.images[n])

@static if VERSION < v"1.11"
    function Base.copy(p::Perm)
        imgs = copy(p.images)
        q = typeof(p)(imgs; check = false)
        if isdefined(p, :inv, :sequentially_consistent)
            inv_imgs = copy(@atomic(p.inv).images)
            q⁻¹ = typeof(p)(inv_imgs; check = false)
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
                σ⁻¹ = typeof(σ)(invperm(σ.images); check = false)
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
else
    function Base.copy(p::Perm)
        imgs = copy(p.images)
        q = typeof(p)(imgs; check = false)
        if isdefined(p, :inv, :acquire)
            inv_imgs = copy(p.inv.images)
            q⁻¹ = typeof(p)(inv_imgs; check = false)
            @atomic :release q⁻¹.inv = q
            @atomiconce :release :acquire q.inv = q⁻¹
        end
        return q
    end

    function Base.inv(σ::Perm)
        if !isdefined(σ, :inv, :acquire)
            if isone(σ)
                @atomiconce :release :acquire σ.inv = σ
            else
                σ⁻¹ = typeof(σ)(invperm(σ.images); check = false)
                # this order is important:
                # fuly initialize the "local" inverse first and only then
                # update σ to make the local inverse visible globally
                @atomic :release σ⁻¹.inv = σ
                @atomiconce :release :acquire σ.inv = σ⁻¹
            end
        end
        return σ.inv
    end

    function AP.cycles(σ::Perm)
        if !isdefined(σ, :cycles, :acquire)
            cdec = AP.CycleDecomposition(σ)
            # we can afford producing more than one cycle decomposition
            @atomiconce :release :acquire σ.cycles = cdec
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
