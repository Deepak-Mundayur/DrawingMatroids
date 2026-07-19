module Maximat

using Base.Threads

include("hypergraph.jl")
include("datamat.jl")
include("comparison.jl")

using .Hypergraph: is_inside, inter_size, diff_size, rem_sub, replace,
                   replace_vertices, identify, remove, remove_element, inf_subs, supp
using .DataMat: cyclic_to_hyper, hyper_to_cyclic, printmat
using .Comparison: test_T3, sup_hyper, poset_maxs_part,
                   serial_max_merge, parallel_max_merge,
                   serial_max_poset, parallel_max_poset

export is_inside, inter_size, diff_size, rem_sub, replace, replace_vertices,
       identify, remove, remove_element, inf_subs, supp,
       cyclic_to_hyper, hyper_to_cyclic, printmat,
       test_T3, sup_hyper, poset_maxs_part,
       serial_max_merge, parallel_max_merge,
       serial_max_poset, parallel_max_poset,
       sup_cyclic, solve_int, detect_mat_case, detect_mat,
       comp_leaves, redund_edge, add_edge, maximal_degenerations, show_dep

"""
    sup_cyclic(CF1, CF2, d)

Test whether the matroid represented by `CF1` is greater than or equal to the
one represented by `CF2`.
"""
sup_cyclic(CF1, CF2, d::Integer) =
    sup_hyper(cyclic_to_hyper(CF1, d), cyclic_to_hyper(CF2, d))


_edge_signature(edge) = Tuple(sort!(Int.(collect(edge))))

function _XT_signature(XT::AbstractVector)
    return (
        Set(_edge_signature(edge) for edge in XT[1]),
        Set(_edge_signature(edge) for edge in XT[2]),
    )
end

function solve_int(XT::AbstractVector, r::Integer)
    flag = false
    T2, T3 = copy(XT[1]), copy(XT[2])

    # Case 4.
    for i in eachindex(T2), j in (i + 1):length(T2)
        if length(intersect(T2[i], T2[j])) == 1 &&
           !is_inside(union(T2[i], T2[j]), T3)
            flag = true
            push!(T3, union(T2[i], T2[j]))
        end
    end
    T2, T3 = rem_sub(Any[T2, T3])

    if r >= 3
        # Case 2. Guarded indices make the original in-place logic safe in Julia.
        i = 1
        while i <= length(T3)
            j = 1
            while j <= length(T2)
                if length(intersect(T3[i], T2[j])) > 1 && !issubset(T2[j], T3[i])
                    flag = true
                    merged = union(T3[i], T2[j])
                    T3 = [T3[k] for k in eachindex(T3) if k != i]
                    push!(T3, merged)
                    i = min(i, length(T3))
                end
                j += 1
            end
            i += 1
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    if r == 3
        i = 1
        while i <= length(T2)
            j = i + 1
            while j <= length(T2)
                if length(intersect(T2[i], T2[j])) > 1
                    flag = true
                    merged = union(T2[i], T2[j])
                    T2 = [T2[k] for k in eachindex(T2) if k != i && k != j]
                    push!(T2, merged)
                    j = i + 1
                else
                    j += 1
                end
            end
            i += 1
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    if r >= 4
        i = 1
        while i <= length(T3)
            j = i + 1
            while j <= length(T3)
                intersection = intersect(T3[i], T3[j])
                if length(intersection) > 2 && !is_inside(intersection, T2)
                    flag = true
                    merged = union(T3[i], T3[j])
                    T3 = [T3[k] for k in eachindex(T3) if k != i && k != j]
                    push!(T3, merged)
                    j = i + 1
                else
                    j += 1
                end
            end
            i += 1
        end
        T2, T3 = rem_sub(Any[T2, T3])
    end

    result = Any[T2, T3]

    # `flag` above records attempted repairs. Some attempted repairs are
    # immediately removed by `rem_sub` as redundant, leaving the reduced
    # hypergraph unchanged. Returning `true` in that situation makes the
    # preprocessing fixed-point loop in `comp_leaves` run forever.
    changed = _XT_signature(result) != _XT_signature(XT)
    return changed, result
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

function comp_leaves(HT::AbstractVector, r::Integer; pproc::Bool=false)
    XT, partition = HT
    XT = Any[copy(XT[1]), copy(XT[2])]
    candidates = Any[]

    XT = rem_sub(XT)
    if pproc
        changed = true
        while changed
            changed, XT = solve_int(XT, r)
        end
    end

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
                    comp_leaves(Any[Any[T2, T3bis], partition], r; pproc=pproc))
        end
        if r <= 3
            T2bis = copy(T2)
            push!(T2bis, intersect(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2bis, T3], partition], r; pproc=pproc))
        end
    elseif case == 2
        e1, e2 = T3[i], T2[j]
        if r <= 4
            T3bis = [T3[k] for k in eachindex(T3) if k != i]
            push!(T3bis, union(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2, T3bis], partition], r; pproc=pproc))
        end
        if r <= 2
            identified = identify(HT, intersect(e1, e2))
            append!(candidates, comp_leaves(identified, r; pproc=pproc))
        end
    elseif case == 3
        e1, e2 = T2[i], T2[j]
        if r <= 3
            T2bis = [T2[k] for k in eachindex(T2) if k != i && k != j]
            push!(T2bis, union(e1, e2))
            append!(candidates,
                    comp_leaves(Any[Any[T2bis, T3], partition], r; pproc=pproc))
        end
        if r <= 2
            identified = identify(HT, intersect(e1, e2))
            append!(candidates, comp_leaves(identified, r; pproc=pproc))
        end
    elseif case == 4
        e1, e2 = T2[i], T2[j]
        T3bis = copy(T3)
        push!(T3bis, union(e1, e2))
        append!(candidates,
                comp_leaves(Any[Any[T2, T3bis], partition], r; pproc=pproc))
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
        return remove_element(HT, minimum(edge))
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

function _combinations(n::Integer, k::Integer)
    k < 0 && return Vector{Vector{Int}}()
    k == 0 && return [Int[]]
    k > n && return Vector{Vector{Int}}()

    result = Vector{Vector{Int}}()
    current = Vector{Int}(undef, k)

    function visit(start::Int, depth::Int)
        if depth > k
            push!(result, copy(current))
            return
        end
        last = n - (k - depth)
        for value in start:last
            current[depth] = value
            visit(value + 1, depth + 1)
        end
    end

    visit(1, 1)
    return result
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

show_dep(candidate, v::Integer) = begin
    v > 1 && println("Add ", candidate[2])
    candidate[1]
end

function _process_candidate(candidate, rank::Integer, preprocess::Bool,
                            verbosity::Integer, existing_maxima)
    HT = show_dep(candidate, verbosity - 1)
    leaves = comp_leaves(HT, rank; pproc=preprocess)
    return poset_maxs_part(leaves, existing_maxima, sup_hyper)
end

"""
    maximal_degenerations(CF, d; S=[1,2,3,4], v=0,
                          preprocess=false, n_procs=1)

Compute all maximal matroid degenerations of a rank-at-most-four matroid given
by cyclic flats grouped by rank.

`n_procs` is interpreted as the maximum number of Julia worker threads. Start
Julia with `julia -t N` to make `N` threads available.
"""
function maximal_degenerations(CF::AbstractVector, d::Integer;
                               S::AbstractVector{<:Integer}=[1, 2, 3, 4],
                               v::Integer=0,
                               preprocess::Bool=false,
                               n_procs::Integer=1)
    @assert all(i -> i in 1:4, S) "Wrong input for S"
    @assert length(CF) == 4 && all(level -> level isa AbstractVector &&
                                  all(edge -> edge isa AbstractSet, level), CF) "Wrong input for CF"
    @assert d > 0 && all(e -> 1 <= e <= d,
                         (e for level in CF for edge in level for e in edge)) "Wrong input for d"

    n_procs <= 0 && (n_procs = Threads.nthreads())
    n_procs = max(1, min(n_procs, Threads.nthreads()))

    HT = cyclic_to_hyper(CF, d)
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

        # The partition, rather than 1:d, is the correct current point set when
        # loops or parallel classes are present.
        for indices in _combinations(length(partition), rank)
            edge = Set(minimum(partition[index]) for index in indices)
            if !redund_edge(XT, edge, rank)
                push!(raw_candidates, (add_edge(HT, edge, rank), edge))
            end
        end

        start = time()
        worker = candidate -> _process_candidate(candidate, rank, preprocess,
                                                  v, cumulative_maxima)
        candidate_maxima = _parallel_map(worker, raw_candidates, n_procs)
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
