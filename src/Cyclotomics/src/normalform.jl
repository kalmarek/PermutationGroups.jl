####
#   normalform (reduction to Zumbroich basis)

normalform(α::Cyclotomic) = normalform!(deepcopy(α))

function normalform!(
    α::Cyclotomic{T},
    tmp = Cyclotomic{T,Vector{T}}(conductor(α), coeffs(α));
    basis_forbidden = zumbroich_viacomplement(conductor(α)),
) where {T}
    # @debug "normalizing:" conductor(α) coeffs(α)
    isnormalized(α, BitSet(first(basis_forbidden))) && return α

    copyto!(coeffs(tmp), coeffs(α))
    normalform!(tmp, basis_forbidden = basis_forbidden)
    copyto!(coeffs(α), coeffs(tmp))

    return α
end

function normalform!(
    α::Cyclotomic{T,<:DenseVector};
    basis_forbidden = zumbroich_viacomplement(conductor(α)),
) where {T}

    basis, forbidden = basis_forbidden
    isnormalized(α, BitSet(basis)) && return α
    for fb in forbidden
        for exp in exponents(α)
            exp in basis && continue
            _replace_exponent!(α, exp, fb)
        end
    end
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

function _tmp_for_reduced_embedding(α::Cyclotomic{T}) where T

    all(iszero, exponents(α)) && return Cyclotomic{T, Vector{T}}(1, [α[0]])

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
        @inbounds for (e,v) in α
            tmp[div(e,k)] = α[e]
        end
    else
        @debug "No trivial reduction is possible"
        copyto!(tmp.coeffs, coeffs(α))
    end
    return tmp
end

"""
    reduced_embedding(α::Cyclotomic{T,V}[, m::Integer=1[,
        tmp = _tmp_for_reduced_embedding(normalform!(α))]])
Return the reduced embedding of `α` into `m`-th cyclotomic field.

If temporary element `tmp` is provided, then
 * `α` must already be in its normal form, and
 * `gcd(exponents(α)) == 1`, so that no trivial reduction is possible.
"""
function reduced_embedding(α::Cyclotomic{T,V}, m::Integer=1,
    tmp = _tmp_for_reduced_embedding(normalform!(α))) where {T, V}

    if conductor(tmp) == 1
        res = Cyclotomic{T, V}(conductor(tmp), coeffs(tmp))
        return res
    end

    basis, forbidden = Cyclotomics.zumbroich_viacomplement(conductor(tmp))

    normalform!(tmp, basis_forbidden=(basis, forbidden))

    phi_nc = length(basis)
    nz = count(!iszero, coeffs(tmp))

    if all(p == q for (p,q,_) in forbidden) && # conductor(tmp) is square-free
        phi_nc == nz && # tmp is supported on φ(n) elements
        _all_equal(tmp, exponents(tmp)) # all nz coeffs of tmp are equal

        val = if iseven(length(forbidden))
            tmp[first(exponents(tmp))]
        else
            -tmp[first(exponents(tmp))]
        end
        @debug "Cyclotomic is real:" α = val
        res = Cyclotomic{T, V}(m, similar(coeffs(α), m))
        zero!(res)
        res[0] = val
        return res
    end

    for (p,q,fb_res) in forbidden
        n = conductor(tmp)
        nz = count(!iszero, coeffs(tmp))

        p == q || continue
        @debug "$(p)² does NOT divide conductor $(conductor(tmp))" p q

        m % p != 0 || continue
        @debug "prime $p does not divide $m so embedding is possible"

        # there are n/p residue classes, each contains p-1 elts, and tmp is constant on each
        nnz, _r = divrem(nz, p-1)
        @debug "the number of non-zero coeffs $nz ≡ $_r mod $(p-1)"
        _r == 0 || continue

        # Check that coeffs(tmp) are constant on each congruence classes.
        # To illustrate by example write:
        #
        # 1 = - ζ₅ - ζ₅² - ζ₅³ - ζ₅⁴
        #
        # and consider
        #
        # c*ζ₉^e = c*ζ₉^e * 1 = c*ζ₉^e * (- ζ₅ - ζ₅² - ζ₅³ - ζ₅⁴) =
        # - c*ζ₄₅^(5e+1*9) - c*ζ₄₅^(5e+2*9) - c*ζ₄₅^(5e+3*9) - c*ζ₄₅^(5e+4*9)
        #
        # If tmp ∈ ℚ(ζ_n) can be reduced to ℚ(ζ_{n÷p})
        # (n = 45, p = 5 in the example), we need to check that all coefficients
        #
        # tmp[p*e + i*n÷p] are equal (for 1 ≤ i ≤ p-1)
        #
        # where exponents 0 ≤ e ≤ n-1 represent congruence classes modulo n÷p

        n_p = n ÷ p

        equal_on_classes = true
        for pe in 0:p:(n-1)
            equal_on_classes || break
            if !_all_equal(tmp, pe.+(n_p.*(1:p-1)))
                equal_on_classes = false
                @debug "skipping reduction over prime $p as cyclotomic is not constant over class $(pe + (n_p)) (mod $(n_p)) " (pe + n_p, tmp[pe + n_p]) (pe.+(n_p.*(2:p-1)), tmp[pe.+(n_p.*(2:p-1))])
            end
        end
        equal_on_classes || continue # to the next prime

        @debug "cyclotomic is constant over residue classes (mod $(n_p)), reducing"

        for pe in 0:p:(n-1)
            # replace sum_i(tmp[p*e + i*n÷p]) (1 ≤ i ≤ p-1)
            # inverse to normalform! over p
            @debug "replacing $pe-th power by $(pe + n_p)-th:" tmp[pe] -tmp[pe + n_p]
            tmp[pe] = -tmp[pe + n_p] # tmp[pe + 0*n_p] = -tmp[pe + 1*n_p]
            tmp[pe.+(n_p.*(1:p-1))] .= zero(T)
        end

        for i in 1:n_p
            tmp[i] = tmp[i*p]
            tmp[i*p] = zero(T)
        end

        tmp = Cyclotomic{T, Vector{T}}(n_p, resize!(tmp.coeffs, n_p))
    end

    res = Cyclotomic{T, V}(conductor(tmp), coeffs(tmp))

    return res
end
