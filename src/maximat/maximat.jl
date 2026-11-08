module Maximat

using Base.Threads
using Combinatorics: combinations

include("hypergraph.jl")
include("datamat.jl")
include("comparison.jl")

using .Hypergraph: is_inside, inter_size, rem_sub, identify, remove, supp
using .DataMat: cyclic_to_hyper, hyper_to_cyclic
using .Comparison: sup_hyper, poset_maxs_part, parallel_max_poset

"""
    sup_cyclic(CF1, CF2, d)

Test whether the matroid represented by `CF1` is greater than or equal to the
one represented by `CF2`.
"""
sup_cyclic(CF1, CF2, d::Integer) =
    sup_hyper(cyclic_to_hyper(CF1, d), cyclic_to_hyper(CF2, d))

"""
    solve_int(XT, r)

Fix all submodularity failures of the labeled hypergraph `XT` that require no choice,
depending on the situation given by `r` (no cyclic flat of rank < `r` included). Return
`(flag, [T2, T3])`, where `flag` tells whether modifications have been done.
"""
function solve_int(XT::AbstractVector, r::Integer)
    flag = false
    T2, T3 = copy(XT[1]), copy(XT[2])
    # Fix Case 4
    for i in eachindex(T2), j in (i + 1):length(T2)
        if length(intersect(T2[i], T2[j])) == 1 && !is_inside(union(T2[i], T2[j]), XT[2])
            flag = true
            push!(T3, union(T2[i], T2[j]))
        end
    end
    T2, T3 = rem_sub(Any[T2, T3])

    if r >= 3
        # Fix Case 2
        for i in 1:length(T3), j in 1:length(T2)
            if length(intersect(T3[i], T2[j])) > 1 && !issubset(T2[j], T3[i])
                flag = true
                T3 = [[T3[k] for k in eachindex(T3) if k != i]; [union(T3[i], T2[j])]]
            end
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    if r == 3
        # Fix Case 3
        for i in 1:length(T2), j in (i + 1):length(T2)
            if length(intersect(T2[i], T2[j])) > 1
                flag = true
                T2 = [[T2[k] for k in eachindex(T2) if k != i && k != j]; [union(T2[i], T2[j])]]
            end
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    if r >= 4
        # Fix Case 1
        for i in 1:length(T3), j in (i + 1):length(T3)
            S = intersect(T3[i], T3[j])
            if length(S) > 2 && !is_inside(S, T2)
                flag = true
                T3 = [[T3[k] for k in eachindex(T3) if k != i && k != j]; [union(T3[i], T3[j])]]
            end
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    return flag, Any[T2, T3]
end


function detect_mat_case(XT::AbstractVector, c::Integer, r::Integer)
    T2, T3 = XT

    if c == 1
        for i in eachindex(T3), j in (i + 1):length(T3)
            intersection = intersect(T3[i], T3[j])
            if length(intersection) > 2 && !is_inside(intersection, T2)
                return 1, (i, j)
            end
        end
    elseif c == 2
        for i in eachindex(T3), j in eachindex(T2)
            if length(intersect(T3[i], T2[j])) > 1 && !issubset(T2[j], T3[i])
                return 2, (i, j)
            end
        end
    elseif c == 3 && r <= 3
        for i in eachindex(T2), j in (i + 1):length(T2)
            length(intersect(T2[i], T2[j])) > 1 && return 3, (i, j)
        end
    elseif c == 4 && r <= 3
        for i in eachindex(T2), j in (i + 1):length(T2)
            if length(intersect(T2[i], T2[j])) == 1 &&
               !is_inside(union(T2[i], T2[j]), T3)
                return 4, (i, j)
            end
        end
    end

    return 0, ()
end

function detect_mat(XT::AbstractVector, r::Integer)
    for c in (3, 2, 1, 4)
        case, indices = detect_mat_case(XT, c, r)
        case != 0 && return case, indices
    end
    return 0, ()
end

function comp_leaves(HT::AbstractVector, r::Integer)
    XT, partition = HT
    XT = Any[copy(XT[1]), copy(XT[2])]
    candidates = Any[]

    # We first reduce the hypergraph as we added some edges
    XT = rem_sub(XT)
    # Former (disabled) preprocessing option: call solve_int before any recursive call.
    # This does not change the number of candidates, see the doc of solve_int for further details.
    # flag = true
    # while flag
    #     flag, XT = solve_int(XT, r)
    # end
    case, indices = detect_mat(XT, r)
    if case == 0
        push!(candidates, Any[XT, partition])
        return candidates
    end

    T2, T3 = XT
    i, j = indices

    if case == 1
        e1, e2 = T3[i], T3[j]
        if r <= 4
            T3bis = [T3[k] for k in eachindex(T3) if k != i && k != j]
            push!(T3bis, union(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2, T3bis], partition], r))
        end
        if r <= 3
            T2bis = copy(T2)
            push!(T2bis, intersect(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2bis, T3], partition], r))
        end
    elseif case == 2
        e1, e2 = T3[i], T2[j]
        if r <= 4
            T3bis = [T3[k] for k in eachindex(T3) if k != i]
            push!(T3bis, union(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2, T3bis], partition], r))
        end
        if r <= 2
            identified = identify(HT, intersect(e1, e2))
            append!(candidates, comp_leaves(identified, r))
        end
    elseif case == 3
        e1, e2 = T2[i], T2[j]
        if r <= 3
            T2bis = [T2[k] for k in eachindex(T2) if k != i && k != j]
            push!(T2bis, union(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2bis, T3], partition], r))
        end
        if r <= 2
            identified = identify(HT, intersect(e1, e2))
            append!(candidates, comp_leaves(identified, r))
        end
    elseif case == 4
        e1, e2 = T2[i], T2[j]
        T3bis = copy(T3)
        push!(T3bis, union(e1, e2))
        append!(candidates,
                comp_leaves(Any[Any[T2, T3bis], partition], r))
    end

    return candidates
end

function redund_edge(XT::AbstractVector, edge::AbstractSet, i::Integer)
    i <= 2 && return false
    i == 3 && return is_inside(edge, XT[1])
    i == 4 && return is_inside(edge, XT[2]) || inter_size(edge, XT[1], 3)
    throw(ArgumentError("i must be in 1:4"))
end

function add_edge(HT::AbstractVector, edge::AbstractSet, i::Integer)
    if i == 1
        return remove(HT, minimum(edge))
    elseif i == 2
        return identify(HT, edge)
    elseif i == 3
        T2 = copy(HT[1][1])
        push!(T2, copy(edge))
        return Any[Any[T2, copy(HT[1][2])], HT[2]]
    elseif i == 4
        T3 = copy(HT[1][2])
        push!(T3, copy(edge))
        return Any[Any[copy(HT[1][1]), T3], HT[2]]
    end
    throw(ArgumentError("i must be in 1:4"))
end

function _parallel_map(f, values::AbstractVector, n_workers::Integer)
    n = length(values)
    n == 0 && return Any[]
    n_workers <= 1 && return [f(value) for value in values]

    workers = min(n_workers, n)
    result = Vector{Any}(undef, n)
    @sync for worker in 1:workers
        Threads.@spawn begin
            for index in worker:workers:n
                result[index] = f(values[index])
            end
        end
    end
    return result
end

_hyper_key(HT) = (Set(Set.(HT[1][1])), Set(Set.(HT[1][2])), Set(Set.(HT[2])))

# Keep only the first occurrence of each hypergraph across all candidate groups.
function _remove_duplicates(groups)
    seen = Set{Any}()
    return [[HT for HT in group if !(_hyper_key(HT) in seen) && (push!(seen, _hyper_key(HT)); true)]
            for group in groups]
end

"""
    maximal_degenerations(CF, d; S=[1,2,3,4], v=0, n_procs=1)

Compute all maximal matroid degenerations of a rank-at-most-four matroid given
by cyclic flats grouped by rank.

`n_procs` is interpreted as the maximum number of Julia worker threads. Start
Julia with `julia -t N` to make `N` threads available.
"""
function maximal_degenerations(CF::AbstractVector, d::Integer;
                               S::AbstractVector{<:Integer}=[1, 2, 3, 4],
                               v::Integer=0,
                               n_procs::Integer=1)
    issubset(S, 1:4) || throw(ArgumentError("Wrong input for S: must be a sublist of [1,2,3,4]"))
    HT = cyclic_to_hyper(CF, d) # Also checks the input CF and d

    n_procs <= 0 && (n_procs = Threads.nthreads())
    n_procs = max(1, min(n_procs, Threads.nthreads()))

    XT, partition = HT
    leaf_time = 0.0
    maxima_time = 0.0

    if v > 0
        names = join(("S$(i)(M)" for i in S), ", ")
        suffix = n_procs > 1 ? "es" : ""
        println("Compute maximals in $names ($n_procs process$suffix)")
    end

    cumulative_maxima = Any[]

    for rank in sort(collect(S); rev=true)
        raw_candidates = Any[]

        # Any non-loop element can become a loop, even inside a parallel class;
        # otherwise take a circuit of size `rank`, given by its representatives.
        edges = rank == 1 ? [Set([e]) for e in sort!(collect(supp(partition)))] :
            [Set(minimum(partition[index]) for index in indices)
             for indices in combinations(1:length(partition), rank)]
        for edge in edges
            if !redund_edge(XT, edge, rank)
                v > 2 && println("Add ", edge)
                push!(raw_candidates, add_edge(HT, edge, rank))
            end
        end

        start = time()
        leaves = _remove_duplicates(_parallel_map(HT -> comp_leaves(HT, rank), raw_candidates, n_procs))
        candidate_maxima = _parallel_map(Y -> poset_maxs_part(Y, cumulative_maxima, sup_hyper), leaves, n_procs)
        after_leaves = time()
        leaf_time += after_leaves - start

        v > 1 && print("Max cands: ", round(after_leaves - start; digits=2),
                       "s ; Intermax: ")

        maximal_groups = parallel_max_poset(candidate_maxima, sup_hyper;
                                             n_procs=n_procs)
        current = Any[]
        for group in maximal_groups
            append!(current, group)
        end
        push!(cumulative_maxima, current)

        maxima_time += time() - after_leaves
        if v > 1
            println(round(time() - after_leaves; digits=2), "s ; ",
                    sum(length, cumulative_maxima; init=0), " current maxs")
        end
    end

    if v > 0
        println("Time elapsed ", round(leaf_time + maxima_time; digits=2),
                "s (total) ; ", round(leaf_time; digits=2),
                "s (max cands) ; ", round(maxima_time; digits=2),
                "s (inter maxs)\n")
    end

    result = Any[]
    for group in cumulative_maxima, hypergraph in group
        push!(result, hyper_to_cyclic(hypergraph, d))
    end
    return result
end

end # module Maximat
