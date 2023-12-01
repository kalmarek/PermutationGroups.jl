#=
function power_by_cycles(σ::AbstractPermutation, n::Integer)
    if n < 0
        return inv(σ)^-n
    elseif n == 0
        return one(σ)
    elseif n == 1
        return copy(σ)
    elseif n == 2
        return σ * σ
    elseif n == 3
        return σ * σ * σ
    elseif n == 4
        return σ * σ * σ * σ
    else
        img = Vector{inttype(σ)}(undef, deg)
        @inbounds for cycle in cycles(σ)
            l = length(cycle)
            k = n % l
            for (idx, j) in enumerate(cycle)
                idx += k
                idx = (idx > l ? idx - l : idx)
                img[j] = cycle[idx]
            end
        end
        return typeof(P)(img, false)
    end
end
=#
