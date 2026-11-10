"""
    draw_poset(M::Matroid; seed=nothing, filename="matroid_poset.png", folder=nothing, edges=:routed, ...)
    draw_poset(Z::RealizationSpace; ...)

Draw every component of the nonbasis variety of `M` at its place in RealizationSpaces'
weak-order poset: `M` on top, more degenerate matroids lower down, and an edge for each
cover relation. `draw_poset(M)` decomposes `nonbasis_variety(M)`; pass that decomposed
variety instead to draw again without recomputing the decomposition.

Nodes are labelled as in RealizationSpaces: `C_j` is `components(variety)[j]`, and `M` marks
the input matroid. With `folder`, each drawn component is also saved there as `C_j.png`.
Within each level, nodes are ordered to sit near the nodes they are joined to. With
`edges=:routed`, each cover relation is a straight line, coloured by its upper node, that bends
through the gap between nodes of any level it passes; with `edges=:labels`, each caption
instead lists the nodes directly above it.

Returns a `PosetDrawing` with fields `plot`, `poset`, `variety`, `labels`, `drawings` and
`seed`: `labels[i]` and `drawings[i]` belong to poset node `i`; a drawing is the
`draw_matroid(C)` result, or a string saying why the node was not drawn. Other keywords
(`attempts`, `tol`, `optimizer`, ...) go to `draw_matroid(C)`.
"""
function draw_poset(M::Oscar.Matroid; kwargs...)
    # Without nonbases the nonbasis variety has no equations: the poset is M alone.
    isempty(Oscar.nonbases(M)) && return _draw_poset(M, nothing; kwargs...)
    return draw_poset(RealizationSpaces.nonbasis_variety(M; compute_matroid_orbits=false); kwargs...)
end

function draw_poset(Z::RealizationSpaces.RealizationSpace; kwargs...)
    RealizationSpaces.mode(Z) in (:nonbasis_variety, :circuit_variety) || throw(ArgumentError(
        "draw_poset needs a decomposed RealizationSpaces.nonbasis_variety(M); a realization space " *
        "keeps only the components of M itself."))
    return _draw_poset(RealizationSpaces.matroid(Z), Z; kwargs...)
end

function _draw_poset(M, Z; seed=nothing, filename="matroid_poset.png", folder=nothing,
                     edges::Symbol=:routed, show_labels::Bool=true, kwargs...)
    edges in (:routed, :labels) || throw(ArgumentError("edges must be :routed or :labels."))
    seed = Int(something(seed, rand(Random.RandomDevice(), UInt32)))
    comps = isnothing(Z) ? RealizationSpaces.RealizationComponent[] : RealizationSpaces.components(Z)

    # Each component's stored matroid, checked before the poset is built from them, or why the
    # check failed.
    checked = map(comps) do C
        try
            _component_matroid!(C; tol=1e-7)
        catch err
            err isa InterruptException && rethrow()
            sprint(showerror, err)
        end
    end
    poset = isnothing(Z) ? RealizationSpaces.weak_order_poset(M, Oscar.Matroid[]) : RealizationSpaces.weak_order_poset(Z)

    # Node i comes from candidate j: component C_j, or M itself when appended after them.
    sources = [only(poset.element_indices[i]) for i in eachindex(poset.elements)]
    labels = map(eachindex(sources)) do i
        j = sources[i]
        j > length(comps) ? "M" : j in poset.input_indices ? "M = C_$j" : "C_$j"
    end

    drawings = map(enumerate(sources)) do (i, j)
        println("[draw_poset] Drawing $(labels[i]) ($i/$(length(sources)))")
        try
            if j > length(comps)
                isempty(Oscar.nonbases(M)) || return "no component of the nonbasis variety has this matroid"
                return draw_matroid(M; kwargs..., seed=seed, filename=nothing)
            end
            checked[j] isa AbstractString && return checked[j]
            return _draw_component(comps[j], checked[j]; kwargs..., seed=seed, filename=nothing)
        catch err
            err isa InterruptException && rethrow()
            return sprint(showerror, err)
        end
    end

    if !isnothing(folder)
        mkpath(folder)
        for (j, d) in zip(sources, drawings)
            d isa AbstractString || savefig(d.plot, joinpath(folder, (j > length(comps) ? "M" : "C_$j") * ".png"))
        end
    end
    plt = _render_poset(poset, labels, drawings; edges, show_labels)
    isnothing(filename) || savefig(plt, String(filename))
    return PosetDrawing(plt, poset, Z, labels, drawings, seed)
end

# The result of draw_poset; the REPL shows a short summary.
struct PosetDrawing
    plot::Plot
    poset
    variety
    labels::Vector{String}
    drawings::Vector
    seed::Int
end

Base.show(io::IO, p::PosetDrawing) = print(io, "PosetDrawing(", length(p.labels), " nodes, seed=", p.seed, ")")

function Base.show(io::IO, ::MIME"text/plain", p::PosetDrawing)
    levels = length(unique(RealizationSpaces._weak_order_levels(p.poset.covers)))
    drawn = count(d -> d isa MatroidDrawing, p.drawings)
    println(io, "Poset drawing: $(length(p.labels)) nodes on $levels levels, $drawn drawn, seed=$(p.seed)")
    print(io, "  fields: plot, poset, variety, labels, drawings, seed")
    for (label, d) in zip(p.labels, p.drawings)
        d isa AbstractString && print(io, "\n  $label: not drawn: ", first(first(split(d, '\n')), 120))
    end
end

# Break text into lines of at most `width` characters, between words.
_wrap(s, width) = rstrip(replace(s, Regex("(.{1,$width})(\\s+|\$)") => s"\1\n"))

# Node positions: RealizationSpaces' poset levels as rows, bottom to top. Within a row, nodes
# are sorted by the mean position of their neighbours in the row above, then in the row below
# (barycenter ordering, a few rounds), which removes most edge crossings.
function _poset_layout(poset; dx=1.6, level_gap=2.4)
    levels = RealizationSpaces._weak_order_levels(poset.covers)
    rows = [findall(==(l), levels) for l in sort(unique(levels))]
    pairs = RealizationSpaces._cover_pairs(poset)                  # (i, j): node i lies below node j
    above = [[j for (k, j) in pairs if k == i] for i in eachindex(levels)]
    below = [[k for (k, j) in pairs if j == i] for i in eachindex(levels)]
    X = zeros(length(levels))
    place!(row) = (X[row] = dx .* ((1:length(row)) .- (length(row) + 1) / 2))
    center(i, nodes) = isempty(nodes[i]) ? X[i] : sum(X[nodes[i]]) / length(nodes[i])
    foreach(place!, rows)
    for _ in 1:4
        foreach(r -> place!(sort!(rows[r]; by=i -> center(i, above))), length(rows)-1:-1:1)
        foreach(r -> place!(sort!(rows[r]; by=i -> center(i, below))), 2:length(rows))
    end
    Y = zeros(length(levels))
    foreach(((r, row),) -> (Y[row] .= (r - 1) * level_gap), enumerate(rows))
    return X, Y
end

# Cover edges as straight lines from below the upper node's caption to above the lower node's
# loop markers; an edge passing a row bends through the gap between boxes nearest to its line.
# Edges of one upper node share a colour.
function _routed_edges!(plt, poset, X, Y, caption_gap)
    pairs = RealizationSpaces._cover_pairs(poset)
    uppers = sort(unique(last.(pairs)))
    colors = palette(:tab10)
    for (i, j) in pairs
        xs, ys = [X[j]], [Y[j] - 0.5 - caption_gap]
        for y in sort(unique(Y); rev=true)
            Y[i] < y < Y[j] || continue
            row = sort(X[Y .== y])
            gaps = [row[1] - 0.8; (row[1:end-1] .+ row[2:end]) ./ 2; row[end] + 0.8]
            x = X[j] + (X[i] - X[j]) * (Y[j] - y) / (Y[j] - Y[i])
            g = gaps[argmin(abs.(gaps .- x))]
            push!(xs, g, g); push!(ys, y + 0.8, y - 0.5 - caption_gap)
        end
        push!(xs, X[i]); push!(ys, Y[i] + 0.75)
        plot!(plt, xs, ys; color=colors[mod1(findfirst(==(j), uppers), length(colors))], linewidth=1.2, label="")
    end
    return plt
end

# Each node is a unit box with a caption below it.
function _render_poset(poset, labels, drawings; edges::Symbol, show_labels::Bool)
    X, Y = _poset_layout(poset)
    caption_gap = 0.45
    xlims = (minimum(X) - 0.8, maximum(X) + 0.8)
    ylims = (minimum(Y) - 0.6 - caption_gap, maximum(Y) + 0.7)
    plt = plot(; legend=:outerbottom, grid=false, axis=false, ticks=false, aspect_ratio=:equal, xlims, ylims, dpi=200,
               size=(max(500, round(Int, 175 * (xlims[2] - xlims[1]))), max(450, round(Int, 175 * (ylims[2] - ylims[1])))))

    edges == :routed && _routed_edges!(plt, poset, X, Y, caption_gap)

    for (i, Mi) in enumerate(poset.elements)
        is_M = startswith(labels[i], "M")
        bounds = ((X[i] - 0.5, X[i] + 0.5), (Y[i] - 0.5, Y[i] + 0.5))
        frame = is_M ? :red : :gray40
        plot!(plt, Shape(X[i] .+ [-0.5, 0.5, 0.5, -0.5], Y[i] .+ [-0.5, -0.5, 0.5, 0.5]);
              fillcolor=:white, linecolor=frame, linewidth=1, label="")

        nonbases = sort(sort.(Oscar.nonbases(Mi)))
        caption = labels[i] * ": rank $(Oscar.rank(Mi)), $(length(nonbases)) nonbases"
        if edges == :labels
            above = [labels[j] for (k, j) in RealizationSpaces._cover_pairs(poset) if k == i]
            isempty(above) || (caption *= "\nbelow: " * join(above, ", "))
        end
        d = drawings[i]
        if d isa AbstractString
            message = "nonbases: " * join(["{" * join(B, ",") * "}" for B in nonbases], " ") * "\n" * first(d, 160)
            annotate!(plt, X[i], Y[i], text(_wrap(message, 34), 5, :black, :center))
        else
            reduction = drawing_reduction(Mi)
            coords = realization_matrix_to_coordinates(d.matrix, reduction; bounds=bounds, real_tol=1e-7)
            _draw_on!(plt, coords, reduction; bounds, show_boundary=false, show_labels,
                      markersize=5, fontsize=5, legend=false)
        end
        annotate!(plt, X[i], Y[i] - 0.5 - caption_gap / 2, text(caption, 7, :black, :center))
    end

    # One legend for all nodes (empty series carry the labels).
    scatter!(plt, [NaN], [NaN]; markershape=:circle, markersize=5, markercolor=:blue, markerstrokecolor=:blue,
             label="points / parallel classes")
    any(Mi -> !isempty(Oscar.loops(Mi)), poset.elements) &&
        scatter!(plt, [NaN], [NaN]; markershape=:rect, markersize=6, markercolor=:orange, markerstrokecolor=:orange,
                 label="loops")
    return plt
end
