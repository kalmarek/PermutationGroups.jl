function __degree(images::AbstractVector{<:Integer})
    deg = 1
    for i in lastindex(images):-1:firstindex(images)
        if images[i] ≠ i
            deg = i
            break
        end
    end
    return deg
end
