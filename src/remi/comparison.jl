module Comparison

using Base.Threads
using ..Hypergraph: supp, inf_subs, replace_vertices, remove_element,
                    is_inside, diff_size

export test_T3, sup_hyper, poset_maxs_part, serial_max_merge,
       parallel_max_merge, serial_max_poset, parallel_max_poset

function test_T3(edge, X2)
    is_inside(edge, X2[2]) && return true
    diff_size(edge, X2[1], 1) && return true
    return false
end

"""
    sup_hyper(H1, H2)

Test whether the matroid encoded by `H1` is greater than or equal to the
matroid encoded by `H2`. Both inputs must be reduced labeled hypergraphs.
"""
function sup_hyper(H1, H2)
    S1, S2 = supp(H1[2]), supp(H2[2])
    issubset(S2, S1) || return false

    H1S2 = deepcopy(H1)
    for element in setdiff(S1, S2)
        H1S2 = remove_element(H1S2, element)
    end

    X1, P1 = H1S2
    X2, P2 = H2
    inf_subs(P1, P2) || return false

    X1bis = replace_vertices(X1, P2)
    inf_subs(X1bis[1], X2[1]) || return false

    return all(edge -> test_T3(edge, X2), X1bis[2])
end

"""
    poset_maxs_part(Y, LX, f)

Keep precisely those maximal elements of `Y` that are not dominated by any
already-known maximal element in the groups `LX`.
"""
function poset_maxs_part(Y::AbstractVector, LX::AbstractVector, f)
    is_maximal = trues(length(Y))

    for i in eachindex(Y)
        for X in LX
            is_maximal[i] || break
            if any(x -> f(x, Y[i]), X)
                is_maximal[i] = false
            end
        end
    end

    for i in eachindex(Y)
        is_maximal[i] || continue
        for j in (i + 1):length(Y)
            is_maximal[j] || continue
            if f(Y[j], Y[i])
                is_maximal[i] = false
                break
            elseif f(Y[i], Y[j])
                is_maximal[j] = false
            end
        end
    end

    return [Y[i] for i in eachindex(Y) if is_maximal[i]]
end

function _split_halves(L::AbstractVector)
    midpoint = length(L) ÷ 2
    return [L[1:midpoint], L[(midpoint + 1):end]]
end

function serial_max_merge(LB::AbstractVector, LC::AbstractVector, f)
    maximal_B = [trues(length(group)) for group in LB]
    maximal_C = [trues(length(group)) for group in LC]

    for i in eachindex(LB)
        for j in eachindex(LB[i])
            maximal_B[i][j] || continue
            for i1 in eachindex(LC)
                maximal_B[i][j] || break
                for j1 in eachindex(LC[i1])
                    maximal_C[i1][j1] || continue
                    if f(LC[i1][j1], LB[i][j])
                        maximal_B[i][j] = false
                        break
                    elseif f(LB[i][j], LC[i1][j1])
                        maximal_C[i1][j1] = false
                    end
                end
            end
        end
    end

    new_B = [[LB[i][j] for j in eachindex(LB[i]) if maximal_B[i][j]]
             for i in eachindex(LB)]
    new_C = [[LC[i][j] for j in eachindex(LC[i]) if maximal_C[i][j]]
             for i in eachindex(LC)]

    return filter(x -> !isempty(x), new_B), filter(x -> !isempty(x), new_C)
end

_can_spawn(depth::Integer, n_procs::Integer) =
    n_procs > 1 && depth < floor(Int, log2(n_procs)) - 1

function parallel_max_merge(LB::AbstractVector, LC::AbstractVector, f;
                            depth::Integer=0, n_procs::Integer=1)
    threshold = 6
    nB = sum(length, LB; init=0)
    nC = sum(length, LC; init=0)

    # The recursive split is over groups, not individual elements. Splitting
    # a zero- or one-group list does not reduce it (`[[], original]`) and
    # therefore recurses forever when that group contains more than the
    # threshold number of elements. Fall back to the exact serial merge.
    if length(LB) <= 1 || length(LC) <= 1 ||
       (nB <= threshold && nC <= threshold)
        return serial_max_merge(LB, LC, f)
    elseif nB > threshold && nC > threshold
        split_B = _split_halves(LB)
        split_C = _split_halves(LC)

        if _can_spawn(depth, n_procs)
            task1 = Threads.@spawn parallel_max_merge(split_B[1], split_C[1], f;
                                                       depth=depth + 1,
                                                       n_procs=n_procs)
            task2 = Threads.@spawn parallel_max_merge(split_B[2], split_C[2], f;
                                                       depth=depth + 1,
                                                       n_procs=n_procs)
            split_B[1], split_C[1] = fetch(task1)
            split_B[2], split_C[2] = fetch(task2)
        else
            split_B[1], split_C[1] = parallel_max_merge(split_B[1], split_C[1], f;
                                                         depth=depth,
                                                         n_procs=n_procs)
            split_B[2], split_C[2] = parallel_max_merge(split_B[2], split_C[2], f;
                                                         depth=depth,
                                                         n_procs=n_procs)
        end

        if _can_spawn(depth, n_procs)
            task1 = Threads.@spawn parallel_max_merge(split_B[1], split_C[2], f;
                                                       depth=depth + 1,
                                                       n_procs=n_procs)
            task2 = Threads.@spawn parallel_max_merge(split_B[2], split_C[1], f;
                                                       depth=depth + 1,
                                                       n_procs=n_procs)
            split_B[1], split_C[2] = fetch(task1)
            split_B[2], split_C[1] = fetch(task2)
        else
            split_B[1], split_C[2] = parallel_max_merge(split_B[1], split_C[2], f;
                                                         depth=depth,
                                                         n_procs=n_procs)
            split_B[2], split_C[1] = parallel_max_merge(split_B[2], split_C[1], f;
                                                         depth=depth,
                                                         n_procs=n_procs)
        end

        return vcat(split_B[1], split_B[2]), vcat(split_C[1], split_C[2])
    elseif nB > threshold
        split_B = _split_halves(LB)
        split_B[1], LC = parallel_max_merge(split_B[1], LC, f;
                                             depth=depth, n_procs=n_procs)
        split_B[2], LC = parallel_max_merge(split_B[2], LC, f;
                                             depth=depth, n_procs=n_procs)
        return vcat(split_B[1], split_B[2]), LC
    else
        split_C = _split_halves(LC)
        LB, split_C[1] = parallel_max_merge(LB, split_C[1], f;
                                             depth=depth, n_procs=n_procs)
        LB, split_C[2] = parallel_max_merge(LB, split_C[2], f;
                                             depth=depth, n_procs=n_procs)
        return LB, vcat(split_C[1], split_C[2])
    end
end

function serial_max_poset(LX::AbstractVector, f)
    is_maximal = [trues(length(group)) for group in LX]

    for i in eachindex(LX)
        for j in eachindex(LX[i])
            is_maximal[i][j] || continue
            for i1 in (i + 1):length(LX)
                is_maximal[i][j] || break
                for j1 in eachindex(LX[i1])
                    is_maximal[i1][j1] || continue
                    if f(LX[i1][j1], LX[i][j])
                        is_maximal[i][j] = false
                        break
                    elseif f(LX[i][j], LX[i1][j1])
                        is_maximal[i1][j1] = false
                    end
                end
            end
        end
    end

    result = [[LX[i][j] for j in eachindex(LX[i]) if is_maximal[i][j]]
              for i in eachindex(LX)]
    return filter(x -> !isempty(x), result)
end

function parallel_max_poset(LX::AbstractVector, f;
                            depth::Integer=0, n_procs::Integer=1)
    threshold = 6
    total = sum(length, LX; init=0)

    # As in `parallel_max_merge`, recursion must strictly shrink the number
    # of groups. A one-group list with more than `threshold` elements was
    # previously split into `[[], LX]`, causing infinite recursion.
    if length(LX) <= 1 || total <= threshold
        return serial_max_poset(LX, f)
    end

    halves = _split_halves(LX)
    if _can_spawn(depth, n_procs)
        task1 = Threads.@spawn parallel_max_poset(halves[1], f;
                                                  depth=depth + 1,
                                                  n_procs=n_procs)
        task2 = Threads.@spawn parallel_max_poset(halves[2], f;
                                                  depth=depth + 1,
                                                  n_procs=n_procs)
        halves[1], halves[2] = fetch(task1), fetch(task2)
    else
        halves[1] = parallel_max_poset(halves[1], f;
                                       depth=depth + 1, n_procs=n_procs)
        halves[2] = parallel_max_poset(halves[2], f;
                                       depth=depth + 1, n_procs=n_procs)
    end

    left, right = parallel_max_merge(halves[1], halves[2], f;
                                     depth=depth, n_procs=n_procs)
    return vcat(left, right)
end

end # module Comparison
