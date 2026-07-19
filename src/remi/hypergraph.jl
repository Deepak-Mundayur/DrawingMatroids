module Hypergraph

import Base: replace

export is_inside, inter_size, diff_size, rem_sub, replace, replace_vertices,
       identify, remove, remove_element, inf_subs, supp

"""
    is_inside(s, L)

Return `true` when the set `s` is contained in at least one set in `L`.
"""
is_inside(s::AbstractSet, L::AbstractVector) = any(l -> issubset(s, l), L)

"""
    inter_size(s, L, n)

Return `true` when `s` intersects some set in `L` in at least `n` elements.
"""
inter_size(s::AbstractSet, L::AbstractVector, n::Integer) =
    any(l -> length(intersect(s, l)) >= n, L)

"""
    diff_size(s, L, n)

Return `true` when `s : l` has exactly `n` elements for some `l` in `L`.
"""
diff_size(s::AbstractSet, L::AbstractVector, n::Integer) =
    any(l -> length(setdiff(s, l)) == n, L)

"""
    rem_sub(XT)

Remove redundant subsets from a two-level labeled hypergraph.
`XT[1]` stores type-2 edges and `XT[2]` stores type-3 edges.
"""
function rem_sub(XT::AbstractVector)
    @assert length(XT) == 2 "XT must contain exactly two edge lists"

    XT_new = [copy(XT[1]), copy(XT[2])]

    # Criterion 1: remove any edge contained in another edge of the same type.
    for k in eachindex(XT_new)
        old = XT_new[k]
        XT_new[k] = [
            s for (i, s) in pairs(old)
            if !any(j -> j != i && issubset(s, old[j]), eachindex(old))
        ]
    end

    # Criterion 2: remove a type-3 edge that is exactly one element larger
    # than a type-2 edge contained in it.
    XT_new[2] = [
        s for s in XT_new[2]
        if !any(t -> length(t) + 1 == length(s) && issubset(t, s), XT_new[1])
    ]

    return XT_new
end

"""
    replace_vertices(XT, P)

Replace each element by the minimum representative of its block in the
partition `P`. Edges that become too small are discarded.
"""
function replace_vertices(XT::AbstractVector, P::AbstractVector)
    representatives = Dict{Int, Int}()
    for block in P
        isempty(block) && continue
        representative = minimum(block)
        for element in block
            representatives[element] = representative
        end
    end

    result = Vector{Any}(undef, length(XT))
    for i in eachindex(XT)
        threshold = i + 2 # 3 for type-2 edges, 4 for type-3 edges
        edges = Set{Int}[]
        for edge in XT[i]
            replaced = Set(representatives[j] for j in edge)
            length(replaced) >= threshold && push!(edges, replaced)
        end
        result[i] = edges
    end
    return result
end

"""
    identify(HT, elements)

Identify all partition blocks meeting `elements`, then reduce the edge lists
using the new representatives.
"""
function identify(HT::AbstractVector, elements::AbstractSet)
    @assert length(HT) == 2 "HT must have the form [XT, P]"
    XT, P = HT

    P1 = Set{Int}[]
    merged = Set{Int}()
    for block in P
        if !isempty(intersect(elements, block))
            union!(merged, block)
        else
            push!(P1, copy(block))
        end
    end
    !isempty(merged) && push!(P1, merged)

    return Any[replace_vertices(XT, P1), P1]
end

"""
    remove_element(HT, e)

Adjoin `{e}` as a type-0 edge, equivalently removing `e` from the represented
matroid and reducing the labeled hypergraph.
"""
function remove_element(HT::AbstractVector, e::Integer)
    @assert length(HT) == 2 "HT must have the form [XT, P]"
    XT, P = HT

    block_index = findfirst(block -> e in block, P)
    block_index === nothing && return HT

    block = P[block_index]

    if length(block) == 1
        P1 = [copy(P[j]) for j in eachindex(P) if j != block_index]
        XT1 = Any[Set{Int}[], Set{Int}[]]

        for edge_type in 1:2
            for edge in XT[edge_type]
                if e in edge
                    # Original minimum sizes are 3 and 4 respectively.
                    if length(edge) > edge_type + 2
                        push!(XT1[edge_type], setdiff(edge, Set([Int(e)])))
                    end
                else
                    push!(XT1[edge_type], copy(edge))
                end
            end
        end
        return Any[XT1, P1]
    end

    new_block = setdiff(block, Set([Int(e)]))
    P1 = Set{Int}[]
    for j in eachindex(P)
        push!(P1, j == block_index ? new_block : copy(P[j]))
    end

    # If e was not the representative, the edge labels do not change.
    if minimum(block) != e
        return Any[[copy(XT[1]), copy(XT[2])], P1]
    end

    new_representative = minimum(new_block)
    XT1 = Any[Set{Int}[], Set{Int}[]]
    for edge_type in 1:2
        for edge in XT[edge_type]
            if e in edge
                new_edge = setdiff(edge, Set([Int(e)]))
                push!(new_edge, new_representative)
                push!(XT1[edge_type], new_edge)
            else
                push!(XT1[edge_type], copy(edge))
            end
        end
    end

    return Any[XT1, P1]
end

# API-compatible aliases for the original Python function names.
replace(XT::AbstractVector, P::AbstractVector) = replace_vertices(XT, P)
remove(HT::AbstractVector, e::Integer) = remove_element(HT, e)

"""
    inf_subs(P1, P2)

Return `true` exactly when every set in `P1` is contained in some set in `P2`.
"""
inf_subs(P1::AbstractVector, P2::AbstractVector) =
    all(p1 -> any(p2 -> issubset(p1, p2), P2), P1)

"""
    supp(L)

Return the union of all sets in `L`.
"""
function supp(L::AbstractVector)
    result = Set{Int}()
    for set in L
        union!(result, set)
    end
    return result
end

end # module Hypergraph
