using Oscar

include("src/remi/maximat.jl")
using .Maximat

M = non_fano_matroid()
d = length(M)

CF = [Set{Int}[] for _ in 1:4]

for F in Oscar.cyclic_flats(M)
    T = Set(Int.(collect(F)))
    k = Oscar.rank(M, collect(T))
    0 <= k <= 3 && push!(CF[k + 1], T)
end

for tested_ranks in ([1], [2], [3])
    println("\nTesting S = $tested_ranks")
    flush(stdout)

    result = nothing

    elapsed = @elapsed begin
        result = Maximat.maximal_degenerations(
            CF,
            d;
            S=tested_ranks,
            preprocess=true,
            v=2,
            n_procs=Threads.nthreads(),
        )
    end

    println(
        "Finished S = $tested_ranks in ",
        round(elapsed; digits=4),
        " seconds; found ",
        length(result),
        " degenerations.",
    )
end