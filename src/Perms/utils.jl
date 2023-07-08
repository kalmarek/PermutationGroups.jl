function __degree(images::AbstractVector{<:Integer})
    for i in lastindex(images):-1:firstindex(images)
        images[i] ≠ i && return i
    end
    return firstindex(images)
end
