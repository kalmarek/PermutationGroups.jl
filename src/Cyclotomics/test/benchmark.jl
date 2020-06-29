using Cyclotomics
using Test
using Statistics

@testset "zumbroich perf" begin
    let ran = 25000:30001
        @assert Cyclotomics.zumbroich_plain(last(ran)) ==
            Cyclotomics.zumbroich(last(ran)) ==
            Cyclotomics.zumbroich_direct(last(ran))

        plain = [(@timed sum(first, Cyclotomics.zumbroich_plain(i)))[2] for i in ran];
        @info "By plain loop: mean and variance over $ran" μ=round(mean(plain), sigdigits=4) σ=round(Statistics.var(plain), sigdigits=4)
        @btime Cyclotomics.zumbroich_plain($(last(ran)))


        compl = [(@timed sum(first, Cyclotomics.zumbroich(i)))[2] for i in ran];
        @info "Constructing the complement: mean and variance over $ran" μ=round(mean(compl), sigdigits=4) σ=round(Statistics.var(compl), sigdigits=4)
        @btime Cyclotomics.zumbroich($(last(ran)))

        direc = [(@timed length(Cyclotomics.zumbroich_direct(i)))[2] for i in ran];
        @info "GAP style: mean and variance over $ran" μ=round(mean(direc), sigdigits=4) σ=round(Statistics.var(direc), sigdigits=4)
        @btime Cyclotomics.zumbroich_direct($(last(ran)));

    end
end
