function Base.show(io::IO, G::PermGroup)
    print(io, "PermGroup( ")
    join(io, gens(G), ", ")
    return print(io, " )")
end

function Base.show(io::IO, ::MIME"text/plain", G::PermGroup)
    o = isdefined(G, :stabchain) ? " of order $(order(StabilizerChain(G)))" : ""
    ngen = ngens(G)

    ioc = IOContext(io, :compact => true)

    println(
        ioc,
        "Permutation group on ",
        ngen,
        " generator",
        ngen > 1 ? "s" : "",
        o,
        " generated by",
    )
    return Base.print_array(ioc, gens(G))
end

function Base.show(io::IO, g::Permutation)
    ioc = IOContext(io, :compact => true)
    return show(ioc, AP.perm(g))
end

function Base.show(io::IO, m::MIME"text/plain", g::Permutation)
    ioc = IOContext(io, :compact => true)
    return show(ioc, m, AP.perm(g))
end
