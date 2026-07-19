module DataMat

using ..Hypergraph: identify, remove_element, supp

export cyclic_to_hyper, hyper_to_cyclic, printmat

function _valid_cyclic_flats(CF)
    length(CF) == 4 || return false
    return all(level -> level isa AbstractVector &&
                       all(edge -> edge isa AbstractSet, level), CF)
end

"""
    cyclic_to_hyper(CF, d)

Convert cyclic flats grouped by rank into the reduced labeled-hypergraph
representation used by Maximat.
"""
function cyclic_to_hyper(CF::AbstractVector, d::Integer)
    @assert _valid_cyclic_flats(CF) "Wrong input for CF"
    @assert d > 0 && all(e -> 1 <= e <= d,
                         (e for level in CF for edge in level for e in edge)) "Wrong input for d"

    partition = [Set([i]) for i in 1:d]
    HT = Any[Any[copy(CF[3]), copy(CF[4])], partition]

    if !isempty(CF[1])
        for e in CF[1][1]
            HT = remove_element(HT, e)
        end
    end

    for block in CF[2]
        HT = identify(HT, block)
    end

    return HT
end

"""
    hyper_to_cyclic(HT, d)

Convert a reduced labeled hypergraph back to cyclic flats grouped by rank.
"""
function hyper_to_cyclic(HT::AbstractVector, d::Integer)
    @assert length(HT) == 2 && length(HT[1]) == 2 "Wrong input for HT"
    @assert d > 0 "Wrong input for d"

    XT, partition = HT
    @assert all(e -> 1 <= e <= d,
                (e for level in XT for edge in level for e in edge)) "Wrong input for d"
    @assert all(e -> 1 <= e <= d,
                (e for block in partition for e in block)) "Wrong input for d"

    loops = setdiff(Set(1:d), supp(partition))

    CF = Vector{Any}(undef, 4)
    CF[1] = isempty(loops) ? Set{Int}[] : [loops]
    CF[2] = [copy(block) for block in partition if length(block) > 1]
    CF[3] = copy(XT[1])
    CF[4] = copy(XT[2])
    return CF
end

_format_set(s::AbstractSet) = "{" * join(sort!(collect(s)), ",") * "}"

"""
    printmat(CF, pref="")

Print cyclic flats grouped by rank, matching the Python program's display.
"""
function printmat(CF::AbstractVector, pref::AbstractString="")
    @assert _valid_cyclic_flats(CF) "Wrong input for CF"

    pieces = String[]
    !isempty(CF[1]) && push!(pieces, "T0: " * join(_format_set.(CF[1]), " "))

    rank1 = [_format_set(block) for block in CF[2] if length(block) > 1]
    !isempty(rank1) && push!(pieces, "T1: " * join(rank1, " "))

    !isempty(CF[3]) && push!(pieces, "T2: " * join(_format_set.(CF[3]), " "))
    !isempty(CF[4]) && push!(pieces, "T3: " * join(_format_set.(CF[4]), " "))

    println(pref, isempty(pieces) ? "" : " " * join(pieces, " ; "))
    return nothing
end

end # module DataMat
