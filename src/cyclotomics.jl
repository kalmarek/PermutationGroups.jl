module Cyclotomics

using Primes

export zumbroich



###############################################################################
#
#   Zumproich basis implementation

function J(k::Integer, p::Integer)
    if p == 2
        return ifelse(iszero(k), 0:0, 0:1)
    else
        pm1 = convert(Int, p) - 1
        return ifelse(iszero(k), 1:pm1, -pm1>>1:pm1>>1 )
    end
end

function zumbroich_plain(n::Integer, m::Integer=1)

    # following the description of
    # https://www.gap-system.org/Manuals/doc/ref/chap60_mj.html#X7F52BEA0862E06F2

    @assert iszero(last(divrem(n, m)))
    isone(n) && return [0]

    factors_n = factor(n)
    factors_m = factor(m)

    exps = Vector{Int}[]
    sizehint!(exps, 2length(factors_n))

    for (p, ν) in factors_n
        for k in factors_m[p]:ν-1
            push!(exps, [div(n, p^(k+1))*j for j in J(k,p)])
        end
    end

    w = (sum.(Iterators.product(exps...))) .% n
    return sort!(vec(w))
end
###############################################################################
#
# Forbidden residues, i.e. the complement of the Zumbroich basis

function _forbidden_residues(n, p, ν, q=p^ν)
    k = div(n, q)
    if isodd(p)
        a = (p^(ν-1) - 1) >> 1
        return BitSet((k*(-a+q):k:k*(a+q)) .% q)
    else # iseven(p)
        return BitSet((k*(q >> 1):k:k*(q-1)) .% q)
    end
end

function _zumbroich_complement(n::Integer, )

    # following the wonderful documenting comments at the top of
    # https://github.com/gap-system/gap/blob/master/src/cyclotom.c

    factor_n = factor(n)
    forbidden = [(p,p^ν)=>_forbidden_residues(n, p, ν, p^ν) for (p, ν) in factor_n]
    isone(n) && return [0], forbidden

    exps = Vector{typeof(n)}(undef, totient(factor_n))
    # @show forbidden
    count = 0
    for i in 0:n-1
        for ((_,q), forbidden_residues) in forbidden
            if i % q ∈ forbidden_residues
                @goto next_i
            end
        end
        count += 1
        exps[count] = i
        @label next_i
    end
    @assert count == length(exps)
    return exps, forbidden
end

zumbroich_complement(n::Integer) = first(_zumbroich_complement(n))


function zumbroich_direct(n::Integer)
    isone(n) && return [0]
    basis = [0]

    factor_n = factor(n)

    if iseven(n)
        for k in 1:factor_n[2]-1
            basis = union!(basis, basis .+ (n>>(k+1)))
        end
    end

    for (p, ν) in factor_n
        p == 2 && continue
        for k in 0:ν-1
            ran = div(n, p^(k+1)).*J(k, p)
            new_basis = typeof(basis)()
            sizehint!(new_basis, length(basis)*length(ran))
            for b in basis
                append!(new_basis, ran .+ b)
            end
            basis = new_basis
        end
    end

    for (i, x) in enumerate(basis)
        basis[i] = mod(basis[i], n)
    end
    unique!(basis)
    sort!(basis)
    return basis
end

zumbroich = zumbroich_complement

end # of module Cyclotomics
