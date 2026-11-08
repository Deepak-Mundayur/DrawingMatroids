module DataMat

using ..Hypergraph: remove, identify, supp, is_inside

"""
    check_cyclic_flats(CF, d)

Return `(loops, classes, T2, T3)`, the normalized data of `CF`: the loops, the parallel
classes, and the cyclic flats of rank 2 and 3 written with the representative `minimum(c)`
of each parallel class `c`.

Throw an `ArgumentError` unless `CF[i]` is exactly the list of cyclic flats of rank `i - 1`
of a matroid of rank at most 4 on `1:d`. Each cyclic flat may be given either in full, or in
the reduced form output by `hyper_to_cyclic` (without loops, a parallel class being given by
its minimum).
"""
function check_cyclic_flats(CF, d)
    fail(msg) = throw(ArgumentError("Wrong input cyclic flats: " * msg))

    d isa Integer && !(d isa Bool) && d > 0 || fail("d must be a positive integer")
    CF isa Union{AbstractVector, Tuple} && length(CF) == 4 ||
        fail("CF must be a list of 4 lists (cyclic flats of rank 0, 1, 2, 3)")
    for (i, L) in enumerate(CF)
        L isa Union{AbstractVector, Tuple} && all(c -> c isa AbstractSet, L) ||
            fail("CF[$(i - 1)] must be a list of sets")
        for c in L
            all(e -> e isa Integer && !(e isa Bool) && 1 <= e <= d, c) ||
                fail("$(_format_set(c)) in CF[$(i - 1)] is not a subset of [1..$d]")
        end
    end
    CF = [[Set{Int}(c) for c in L] for L in CF]

    # Rank 0: the set of loops
    length(CF[1]) > 1 && fail("CF[0] must contain at most one set (the loops)")
    loops = isempty(CF[1]) ? Set{Int}() : CF[1][1]

    # Rank 1: the parallel classes (possibly given with the loops)
    classes = Set{Int}[]
    for c in CF[2]
        !isempty(c ∩ loops) && !(loops ⊆ c) && fail("rank-1 flat $(_format_set(c)) contains some loops but not all")
        length(setdiff(c, loops)) < 2 && fail("rank-1 flat $(_format_set(c)) has less than 2 non-loop elements")
        c = setdiff(c, loops)
        k = findfirst(c1 -> !isempty(c ∩ c1), classes)
        k === nothing || fail("rank-1 flats $(_format_set(c)) and $(_format_set(classes[k])) intersect")
        push!(classes, c)
    end

    # Rank 2 and 3: written with the representatives of the parallel classes
    rep = Dict(e => minimum(c) for c in classes for e in c)
    T = Vector{Set{Int}}[]
    for i in (2, 3)
        Ti = Set{Int}[]
        for c in CF[i + 1]
            !isempty(c ∩ loops) && !(loops ⊆ c) && fail("rank-$i flat $(_format_set(c)) contains some loops but not all")
            for p in classes
                cp = c ∩ p
                !isempty(cp) && cp != p && cp != Set([minimum(p)]) &&
                    fail("rank-$i flat $(_format_set(c)) contains only part of the parallel class $(_format_set(p))")
            end
            r = Set(get(rep, e, e) for e in setdiff(c, loops))
            length(r) < i + 1 && fail("rank-$i flat $(_format_set(c)) has less than $(i + 1) elements up to parallelism")
            r in Ti && fail("rank-$i flat $(_format_set(c)) is given twice")
            push!(Ti, r)
        end
        push!(T, Ti)
    end
    T2, T3 = T

    # The flats define a matroid whose rank-2/3 cyclic flats are exactly T2, T3
    # (Lemma 5.2 of the paper, together with the cyclicity of the edges of T3)
    for (i, l1) in enumerate(T2), l2 in T2[(i + 1):end]
        length(l1 ∩ l2) > 1 && fail("rank-2 flats $(_format_set(l1)) and $(_format_set(l2)) share more than one point")
        length(l1 ∩ l2) == 1 && !is_inside(l1 ∪ l2, T3) &&
            fail("intersecting rank-2 flats $(_format_set(l1)) and $(_format_set(l2)) are not in a rank-3 flat")
    end
    for h in T3, l in T2
        length(h ∩ l) > 1 && !(l ⊆ h) && fail("rank-2 flat $(_format_set(l)) meets rank-3 flat $(_format_set(h)) in more than one point")
        length(setdiff(h, l)) <= 1 && fail("rank-3 flat $(_format_set(h)) is not cyclic or has rank <= 2 (see $(_format_set(l)))")
    end
    for (i, h1) in enumerate(T3), h2 in T3[(i + 1):end]
        length(h1 ∩ h2) > 2 && !is_inside(h1 ∩ h2, T2) &&
            fail("rank-3 flats $(_format_set(h1)) and $(_format_set(h2)) meet in a set of rank 3")
    end

    return loops, classes, T2, T3
end

"""
    cyclic_to_hyper(CF, d)

Reduced labeled hypergraph of the cyclic flats `CF` (grouped by rank) on ground set `1:d`.
"""
function cyclic_to_hyper(CF, d)
    loops, classes, T2, T3 = check_cyclic_flats(CF, d)

    HT = Any[Any[T2, T3], [Set([i]) for i in 1:d]]
    for e in sort!(collect(loops))
        HT = remove(HT, e) # Remove loops
    end
    for L in classes
        HT = identify(HT, L) # Identify double points
    end
    return HT
end

"""
    hyper_to_cyclic(HT, d)

Cyclic flats grouped by rank encoded in the reduced hypergraph `HT` with ground set `1:d`.
"""
function hyper_to_cyclic(HT::AbstractVector, d::Integer)
    @assert length(HT) == 2 && length(HT[1]) == 2 &&
            all(t -> t isa AbstractVector && all(c -> c isa AbstractSet, t), HT[1]) &&
            all(t -> t isa AbstractSet, HT[2]) "Wrong input for HT"
    @assert d > 0 && all(s -> s <= d, (s for cf in HT[1] for c in cf for s in c)) &&
            all(c -> c <= d, (c for cf in HT[2] for c in cf)) "Wrong input for d"

    XT, P = HT
    loops = setdiff(Set(1:d), supp(P))
    return Any[isempty(loops) ? Set{Int}[] : [loops],      # Type 0
               [copy(p) for p in P if length(p) > 1],      # Type 1
               XT[1],                                      # Type 2
               XT[2]]                                      # Type 3
end

# Python's display of a set of small integers, e.g. {1, 2, 3}.
_format_set(s) = "{" * join(sort!(collect(s)), ", ") * "}"

"""
    printmat(CF, pref="")

Print the matroid defined by the cyclic flats `CF`: after "Ti:" are listed all rank-i cyclic
flats.
"""
function printmat(CF::AbstractVector, pref::AbstractString="")
    @assert length(CF) == 4 && all(t -> t isa AbstractVector && all(c -> c isa AbstractSet, t), CF) "Wrong input for CF"

    Ci = String[]
    !isempty(CF[1]) && push!(Ci, "T0: {" * join(_format_set.(CF[1]), ",") * "}")
    push!(Ci, "T1: " * join(("{" * join(sort!(collect(p)), ",") * "}" for p in CF[2] if length(p) > 1), " "))
    push!(Ci, "T2: " * join((replace(_format_set(l), " " => "") for l in CF[3]), " "))
    push!(Ci, "T3: " * join((replace(_format_set(l), " " => "") for l in CF[4]), " "))
    filter!(c -> length(c) > 4, Ci)
    println(pref * " " * join(Ci, " ; "))
    return nothing
end

end # module DataMat
