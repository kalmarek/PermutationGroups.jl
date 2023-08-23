import PrecompileTools

PrecompileTools.@setup_workload begin
    S = [perm"(1,3,5,7)(2,4,6,8)", perm"(1,3,8)(4,5,7)"]

    PrecompileTools.@compile_workload begin
        G = PermGroup(S)
        @assert order(G) == order(Int, G) == 24
        ff(g) = 1^g
        @assert sum(ff, G) == 108
    end
end
