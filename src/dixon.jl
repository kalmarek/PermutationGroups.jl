import AbstractAlgebra: exponent

exponent(G::AbstractAlgebra.Group) = exponent(conjugacy_classes(G))
exponent(cclasses::AbstractVector) = lcm(order.(first.(cclasses)))

dixon_prime(G::AbstractAlgebra.Group) = dixon_prime(order(G), exponent(G))

function dixon_prime(cclasses::AbstractVector)
    ordG = sum(length, cclasses)
    m = exponent(cclasses)
    return dixon_prime(ordG, m)
end

function dixon_prime(ordG::Integer, exponent::Integer)
    p = 2floor(Int, sqrt(ordG))
    while true
        p = nextprime(p+1)
        isone(p % exponent) && break # we need -1 to be in the field
    end
    return p
end

Base.size(M::CCMatrix) = size(M.m)
Base.IndexStyle(::Type{<:CCMatrix}) = IndexCartesian()

function Base.getindex(M::CCMatrix, s::Integer, t::Integer)
    if isone(-M.m[s,t])
        M.m[s,:] .= 0
        r = M.r
        out = one(first(first(M.cc)))

        for g in M.cc[r]
            for h in M.cc[s]
                out = mul!(out, g, h)
                for t in 1:size(M, 1)
                    if out == first(M.cc[t])
                        M.m[s, t] += 1
                        break
                    end
                end
            end
        end

    end
    return M.m[s,t]
end
