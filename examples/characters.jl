include("../src/characters.jl")

# Build group ring for the permutation group over 4 elements.
N = 4
G = PermGroup(N)
RG = GroupRing(G)

# Conjugacy classes as all-one elements of GroupRing
# WTH does all-one elements mean?
ccG = let ptypes = [p.part for p in AllParts(N)]
    perms_typed = Dict(g => permtype(g) for g in G)
    ccG = [Set(g for (g,t) in perms_typed if t==pt) for pt in ptypes]
    [sum(RG.(cc)) for cc in ccG]
end

# Multiplication tables for conjugacy classes.
M = [CCMatrix(ccG, i) for i in 1:length(ccG)];

# First class is the conjugacy class of e. Multiplication with e should not change any conjugacy class.
@assert M[1] == one(M[1])

# In the sequal we need to do exact arithmetics which is why we work over a finite field. 
# There are probably smarter ways to find a suitable prime...
m = Int(lcm(order.(G)))
q = find_prime(m)
F = AbstractAlgebra.GF(q)


# hopefully don't need the following
@ncpolyvar T
p = mod(characteristic_polynomial(M[2], T), q)


