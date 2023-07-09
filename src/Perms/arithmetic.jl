function Base.inv(σ::AbstractPermutation)
    img = Vector{inttype(σ)}(undef, degree(σ))
    for i in Base.OneTo(degree(σ))
        @inbounds img[i^σ] = i
    end
    return typeof(σ)(img, false)
end

function Base.:(*)(σ::AbstractPermutation, τ::AbstractPermutation)
    deg = max(degree(σ), degree(τ))
    img = Vector{inttype(σ)}(undef, deg)
    for i in Base.OneTo(deg)
        img[i] = (i^σ)^τ
    end
    return typeof(σ)(img, false)
end

function Base.:(*)(
    σ::AbstractPermutation,
    τ::AbstractPermutation,
    ρ::AbstractPermutation,
)
    deg = max(degree(σ), degree(τ), degree(ρ))
    img = Vector{inttype(σ)}(undef, deg)
    for i in Base.OneTo(deg)
        img[i] = ((i^σ)^τ)^ρ
    end
    return typeof(σ)(img, false)
end

function Base.:(*)(σ::AbstractPermutation, τs::AbstractPermutation...)
    deg = max(degree(σ), maximum(degree, τs))
    img = Vector{inttype(σ)}(undef, deg)
    for i in Base.OneTo(deg)
        j = (i^σ)
        for τ in τs
            j = j^τ
        end
        img[i] = j
    end
    return typeof(σ)(img, false)
end

Base.:^(σ::AbstractPermutation, τ::AbstractPermutation) = conj(σ, τ)

function Base.conj(σ::AbstractPermutation, τ::AbstractPermutation)
    deg = max(degree(σ), degree(τ))
    img = Vector{inttype(σ)}(undef, deg)
    for i in Base.OneTo(deg)
        img[i^τ] = (i^σ)^τ
    end
    P = typeof(σ)
    return P(img, false)
end
