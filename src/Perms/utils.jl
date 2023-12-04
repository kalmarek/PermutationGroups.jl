# function __degree_simple(images::AbstractVector{<:Integer})
#     @inbounds for i in lastindex(images):-1:firstindex(images)
#         images[i] ≠ i && return i
#     end
#     return firstindex(images) - 1
# end

@inline function __unsafe_findlast(v, lastidx, firstidx)
    @inbounds for idx in lastidx:-1:firstidx
        v[idx] ≠ idx && return idx
    end
    return lastidx + 1 # so that we're not returning nothing
end

@inline function __unsafe_findlast_simd(v, lastidx, firstidx)
    found = false
    @inbounds for idx in firstidx:lastidx
        found |= v[idx] ≠ idx
    end
    if found
        return @inbounds __unsafe_findlast(v, lastidx, firstidx)
    else
        return lastidx + 1
    end
end

function __degree(v::AbstractVector{<:Integer})
    step = 64
    i = length(v)
    @inbounds while i > step
        ans = __unsafe_findlast_simd(v, i, i + 1 - step)
        ans ≤ i && return ans
        i -= step
    end
    if i > firstindex(v)
        ans = __unsafe_findlast_simd(v, i, firstindex(v))
        return ifelse(ans ≤ i, ans, firstindex(v) - 1)
    else
        return firstindex(v) - 1
    end
end
