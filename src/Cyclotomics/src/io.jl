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
            coeff = α[exp]
            if iszero(exp)
                print(io, coeff)
                continue
            end
            if isone(coeff)
                sign_str = i == 1 ? "" : "+"
                coeff_str = ""
            elseif isone(-coeff)
                sign_str = "-"
                coeff_str = ""
            elseif coeff > zero(coeff)
                sign_str = "+"
                coeff_str = "$coeff*"
            else
                sign_str = ""
                coeff_str = "$coeff*"
            end

            exp_str = isone(exp) ? "" : "$(superscriptify(exp))"
            print(io, " ", sign_str, coeff_str, ζ, exp_str)
        end
    end
end

function Base.print(io::IO, α::Cyclotomic)
    E_str = "E($(conductor(α)))"
    for exp in exponents(α)
        coeff = α[exp]
        if coeff > zero(coeff)
            sign_str = "+"
        else
            sign_str = ""
        end

        print(io, " ", sign_str, coeff, "*", E_str, "^", exp)
    end
end
