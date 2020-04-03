export square_free_factorization

mutable struct SFF{T}
    facts::Vector{T}
    mults::Vector{Int}
end

SFF{T}() where T= SFF(T[], Int[])

function Base.:^(sff::SFF, e::Int)
    return SFF(sff.facts, sff.mults.*e)
end

function Base.show(io::IO, sff::SFF)
    for i = 1:length(sff.mults)-1
        print(io, "( $(sff.facts[i]) )^$(sff.mults[i]) * ")
    end
    print(io, "( $(sff.facts[end]) )^$(sff.mults[end])")
end

function pth_root(a, T, char::Int)
    @assert derivative(a) == 0
    k = Int((length(a.coeffs)-1)/char)
    b = a.coeffs[1]
    for i = 1:k
        b += a.coeffs[i*char+1]*T^i
    end
    return b
end

function square_free_factorization(a, T, char) 
    @assert last(a.coeffs) == 1 "Input not monic."
    i = 1
    factors = SFF{typeof(a)}()
    b = derivative(a)
    if !( b == 0)
        c = gcd(a, b)

        _, w = divides(a, c)
        while !(w == 1)
            y = gcd(w, c)
            _, z = divides(w, y)
            if !(z == 1)
                push!(factors.facts, z)
                push!(factors.mults, i)
            end
            i += 1
            w = y
            _, c = divides(c, y)
        end
        if !( c == 1)
            c = pth_root(c, T, char)
            nsff = square_free_factorization(c, T, char)^char
            append!(factors.facts, nsff.facts)
            append!(factors.mults, nsff.mults)
        end
    else
        a = pth_root(a, T, char)
        factors = square_free_factorization(a, T, char)^char
    end
    return factors
end
