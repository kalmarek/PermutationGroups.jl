####
#   IO

function subscriptify(n::Int)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join((Char(subscript_0 + i) for i in reverse(digits(n))))
end

function superscriptify(n::Int)
    superscripts = Dict(
    0 => "⁰",
    1 => "¹",
    2 => "²",
    3 => "³",
    4 => "⁴",
    5 => "⁵",
    6 => "⁶",
    7 => "⁷",
    8 => "⁸",
    9 => "⁹",
)
    return join(superscripts[d] for d in reverse(digits(n)))
end

function Base.show(io::IO, α::Cyclotomic{T}) where T
    ζ = "ζ"*subscriptify(conductor(α))
    if iszero(α)
        print(io, zero(T))
    else
        for (i, exp) in enumerate(exponents(α))
            val = α[exp]
            if iszero(exp)
                print(io, val)
                continue
            end
            sign_str = (val >= zero(T) ? (i == 1 ? "" : " +") : " ")
            val_str = "$val"
            exp_str = isone(exp) ? "" : "$(superscriptify(exp))"
            print(io, sign_str, val_str, "*", ζ, exp_str)
        end
    end
end
