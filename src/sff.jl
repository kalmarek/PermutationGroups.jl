"""
   Algorithm: SFF (Square-Free Factorization)
   Input: A monic polynomial f in Fq[x] where q=p^m
   Output: Square-free factorization of f
   R ← 1
   # Make w be the product (without multiplicity) of all factors of f that have 
   # multiplicity not divisible by p
   c ← gcd(f, f′)
   w ← f/c 
   
   # Step 1: Identify all factors in w
   i←1 
   while w ≠ 1 do
       y ← gcd(w, c)
       fac ← w/y
       R ← R·faci
       w ← y; c ← c/y; i ← i+1 
   end while
   # c is now the product (with multiplicity) of the remaining factors of f
   
   # Step 2: Identify all remaining factors using recursion
   # Note that these are the factors of f that have multiplicity divisible by p
   if c ≠ 1 then
       c ← c1/p
       R ← R·SFF(c)p
   end if 
   
   Output(R)
"""
#=
F = GF(2)
M = matrix(F, [1 0; 1 1])
R, T = PolynomialRing(F, string(gensym()))
p = minpoly(R, M)
=#
export square_free_factorization
function square_free_factorization(f, char)    
    factors = typeof(f)[]
    c = gcd(f, derivative(f))
    _, w = divides(f, c)
    while !(w == 1)
        y = gcd(w, c)
        _, fac = divides(w, y)
        push!(factors, fac)
        w = y
        _, c = divides(c, y)
    end
    
    # do we have to care about the second step?
    return factors
end
