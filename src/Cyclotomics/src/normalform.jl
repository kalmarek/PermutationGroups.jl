####
#   normalform (reduction to Zumbroich basis)

normalform(α::Cyclotomic) = normalform!(deepcopy(α))

function normalform!(α::Cyclotomic{T};
    tmp=Cyclotomic{T, Vector{T}}(conductor(α), coeffs(α)),
    basis_forbidden = zumbroich_viacomplement(conductor(α))
    ) where T
    # @debug "normalizing:" conductor(α) coeffs(α)

    basis, forbidden = basis_forbidden

    isnormalized(α, BitSet(basis)) && return α

    tmp.coeffs .= coeffs(α)

    for fb in forbidden
        for exp in exponents(tmp)
            exp in basis && continue
            _replace_exponent!(tmp, exp, fb)
        end
    end

    α.coeffs .= coeffs(tmp)
    return α
end

function _replace_exponent!(α::Cyclotomic, exp::Integer, fb)
    val = α[exp]
    p, q, forbidden_exps = fb

    # @debug p forbidden_exps
    if exp % q ∈ forbidden_exps
        # @debug "removing $(α[exp])*ζ^$exp from" α
        m = conductor(α) ÷ p
        # @debug " " collect(α[exp .+ 1:m:m*(p-1)])

        # TODO: this doesn't work
        # α[exp .+ 1:m:m*(p-1)] .-= val

        for i in 1:p-1
            α[exp + i*m] -= val
            # exp + i*m is not guaranteed to be allowed, but it's larger than exp and hence will be dealt later
        end
        α[exp] = 0
        # @debug "after removal:" α
    end
    return α
end

function embed(α::Cyclotomic, m::Integer)
    conductor(α) == m && return deepcopy(α)
    if conductor(α) > m
        @assert conductor(α) % m == 0 "Embeding of ℚ(ζ$(subscriptify(m))) ↪ ℚ(ζ$(subscriptify(conductor(α)))) is not possible."
        return reduced_embedding(α, m)
    else
        k, _r = divrem(m, conductor(α))
        @assert _r == 0 "Embeding of ℚ(ζ$(subscriptify(conductor(α)))) ↪ ℚ(ζ$(subscriptify(m))) is not possible."

        res = zero!(similar(α, valtype(α), m))

        @inbounds for e in exponents(α)
            res[e*k] = α[e]
        end
        return res
    end
end

function reduced_embedding(α::Cyclotomic{T,V}, m::Integer=1) where {T, V}
    k = gcd(exponents(α)...)
    @debug "gcd(exponents(α)) = $k" collect(exponents(α))

    tmp = let (d,r) = divrem(conductor(α), k)
        if (k > 1 && iszero(r))
            @debug "Reducing the embedding ring to $T(ζ$(subscriptify(d)))"
            tmp = Cyclotomic{T, Vector{T}}(d, zeros(T, d))
        else
            n = conductor(α)
            tmp = Cyclotomic{T, Vector{T}}(n, zeros(T, n))
        end
        tmp # tmp has dense storage
    end

    if k > 1
        # the trivial reduction E(n)^e → E(n÷k)^(e÷k)
        @debug "Performing trivial reduction from ℚ(ζ$(subscriptify(conductor(α)))) → ℚ(ζ$(subscriptify(conductor(tmp))))"
        @inbounds for (e,v) in zip(exponents(α), values(α))
            tmp[div(e,k)] = α[e]
        end
    else
        @debug "No trivial reduction is possible"
        tmp.coeffs .= coeffs(α)
    end

    normalform!(tmp)

    basis, forbidden = Cyclotomics.zumbroich_viacomplement(conductor(tmp))

    phi_nc = length(basis)
    nz = count(!iszero, coeffs(tmp))
    val = first(values(tmp))

    if all(p == q for (p,q,_) in forbidden) && # conductor(tmp) is square-free
        phi_nc == nz && # tmp is supported on φ(n) elements
        _all_equal(tmp, exponents(tmp), val) # all nz coeffs of tmp are equal

        @debug "Cyclotomic is real:" α = -val
        res = Cyclotomic{T, V}(m, similar(coeffs(α), m))
        zero!(res)
        res[0] = -val
        return res
    end

    n = conductor(tmp)
    for (p,q,fb_res) in forbidden

        p != q && continue
        @debug "$(p)² does NOT divide conductor $(conductor(tmp))" p q

        m % p == 0 && continue
        @debug "prime $p does not divide $m so embedding is possible"

        # there are nc/p residue classes, each contains p-1 elts, and tmp is constant on each
        nnz, _r = divrem(nz, p-1)
        @debug "the number of non-zero coeffs $nz ≡ $_r mod $(p-1)"
        _r != 0 && continue


        n╱p = n ÷ p

        # check that coeffs(tmp) are constant on each congruence classes
        equal_on_classes = true
        for i in 1:p:n
            if !_all_equal(tmp, i.+(2n╱p:n╱p:n+i), tmp[i + n╱p])
                equal_on_classes = false
                @debug "skipping reduction over prime $p as cyclotomic is not constant over class $(i+(n╱p)) (mod $(n╱p)) " tmp[i.+(2n╱p:n╱p:n+i)]
                break
            end
        end
        equal_on_classes || continue # to next prime

        @debug "cyclotomic is constant over residue classes (mod $(n╱p)), reducing"


        for i in 1:p:n # inverse to normalform over p
            tmp[i] = -tmp[i+n╱p]
            tmp[(i+n╱p):n╱p:(n╱p+i)] = zero(T)
        end

        for i in 1:n╱p
            tmp[i] = tmp[i*p]
            tmp[i*p] = zero(T)
        end

        nz = nnz # update the number on non-zero coefficients
        tmp = Cyclotomic{T, Vector{T}}(np, resize!(coeffs(tmp), np))
    end


    return Cyclotomic{T, V}(conductor(tmp), coeffs(tmp))
end
