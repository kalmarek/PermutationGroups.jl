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
    while p < 1000
        p = nextprime(p+1)
        iszero((p-1) % exponent) && break
    end
    return p
end

