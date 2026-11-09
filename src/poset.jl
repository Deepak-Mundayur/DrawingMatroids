"""
    draw_poset(M::Matroid; seed=nothing, filename="matroid_poset.png", folder=nothing, per_row=4, edges=:routed, ...)
    draw_poset(Z::RealizationSpace; ...)

Draw every component of the nonbasis variety of `M` at its place in RealizationSpaces'
weak-order poset: `M` on top, more degenerate matroids lower down, and an edge for each
cover relation. `draw_poset(M)` decomposes `nonbasis_variety(M)`; pass that decomposed
variety instead to draw again without recomputing the decomposition.

Nodes are labelled as in RealizationSpaces: `C_j` is `components(variety)[j]`, and `M` marks
the input matroid. With `folder`, each drawn component is also saved there as `C_j.png`.
A level with more than `per_row` nodes is wrapped into several rows. With `edges=:routed`,
each cover relation is a line through the gaps between nodes, coloured by its upper node; with
`edges=:labels`, each caption instead lists the nodes directly above it.

Returns `(plot, poset, variety, labels, drawings, seed)`: `labels[i]` and `drawings[i]` belong
to poset node `i`; a drawing is the `draw_matroid(C)` result, or a string saying why the node
was not drawn. Other keywords (`attempts`, `tol`, optimizer settings) go to `draw_matroid(C)`.
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

function _draw_poset(M, Z; seed=nothing, filename="matroid_poset.png", folder=nothing, per_row::Int=4,
                     edges::Symbol=:routed, show_labels::Bool=true, kwargs...)
    edges in (:routed, :labels) || throw(ArgumentError("edges must be :routed or :labels."))
    seed = Int(something(seed, rand(Random.RandomDevice(), UInt32)))
    components = isnothing(Z) ? RealizationSpaces.RealizationComponent[] : RealizationSpaces.components(Z)

    # Check each component's stored matroid before the poset is built from them.
    unchecked = Dict{Int,String}()
    for (j, C) in enumerate(components)
        try
            _component_matroid!(C; tol=1e-7)
        catch err
            err isa InterruptException && rethrow()
            unchecked[j] = sprint(showerror, err)
        end
    end
    poset = isnothing(Z) ? RealizationSpaces.weak_order_poset(M, Oscar.Matroid[]) : RealizationSpaces.weak_order_poset(Z)

    # Node i comes from candidate j: component C_j, or M itself when appended after them.
    sources = [only(poset.element_indices[i]) for i in eachindex(poset.elements)]
    labels = map(eachindex(sources)) do i
        j = sources[i]
        j > length(components) ? "M" : j in poset.input_indices ? "M = C_$j" : "C_$j"
    end

    drawings = map(sources) do j
        try
            if j > length(components)
                isempty(Oscar.nonbases(M)) || return "no component of the nonbasis variety has this matroid"
                return draw_matroid(M; kwargs..., seed=seed, filename=nothing)
            end
            haskey(unchecked, j) && return unchecked[j]
            return draw_matroid(components[j]; kwargs..., seed=seed, filename=nothing)
        catch err
            err isa InterruptException && rethrow()
            return sprint(showerror, err)
        end
    end

    if !isnothing(folder)
        mkpath(folder)
        for (j, d) in zip(sources, drawings)
            d isa AbstractString || savefig(d.plot, joinpath(folder, (j > length(components) ? "M" : "C_$j") * ".png"))
        end
    end
    plt = _render_poset(poset, labels, drawings; per_row, edges, show_labels)
    isnothing(filename) || savefig(plt, String(filename))
    return (plot=plt, poset=poset, variety=Z, labels=labels, drawings=drawings, seed=seed)
end

# Break text into lines of at most `width` characters, between words.
function _wrap(s, width)
    lines = String[]
    for paragraph in split(s, '\n')
        line = ""
        for word in split(paragraph)
            if !isempty(line) && length(line) + 1 + length(word) > width
                push!(lines, line); line = word
            else
                line = isempty(line) ? word : line * " " * word
            end
        end
        push!(lines, line)
    end
    return join(lines, "\n")
end

# RealizationSpaces' poset levels, bottom to top. A level with more than `per_row` nodes is
# wrapped into stacked sub-rows, packed closer than separate levels; returns the node
# positions and the vertical extent of each wrapped level, which is shaded as one group.
function _poset_layout(covers; per_row=4, dx=1.6, row_gap=2.0, level_gap=2.4)
    levels = RealizationSpaces._weak_order_levels(covers)
    X, Y = zeros(length(levels)), zeros(length(levels))
    wrapped = Tuple{Float64,Float64}[]
    y = 0.0
    for level in sort(unique(levels))
        rows = collect(Iterators.partition(findall(==(level), levels), per_row))
        top = y + (length(rows) - 1) * row_gap
        for (r, row) in enumerate(rows), (k, node) in enumerate(row)
            X[node] = dx * (k - (length(row) + 1) / 2)
            Y[node] = top - (r - 1) * row_gap
        end
        length(rows) > 1 && push!(wrapped, (y, top))
        y = top + level_gap
    end
    return X, Y, wrapped
end

# Cover edges routed through white space only: from below the upper node's caption into the
# gap under its row, across each row in between through the gap between boxes nearest to the
# target, then down into the lower node. Edges of one upper node share a colour and an offset.
function _routed_edges!(plt, poset, X, Y, caption_gap)
    rows = sort(unique(Y); rev=true)
    row_of = Dict(y => r for (r, y) in enumerate(rows))
    bottom(r) = rows[r] - 0.5 - caption_gap
    channel(r) = (bottom(r) + rows[r + 1] + 0.65) / 2            # midway to the next row's loop markers
    gutters = map(rows) do y
        xs = sort(X[Y .== y])
        [xs[1] - 0.8; (xs[1:end-1] .+ xs[2:end]) ./ 2; xs[end] + 0.8]
    end
    pairs = RealizationSpaces._cover_pairs(poset)               # (i, j): node i lies below node j
    uppers = sort(unique(last.(pairs)))
    colors = palette(:tab10)
    for (i, j) in pairs
        k = findfirst(==(j), uppers)
        offset = 0.05 * (mod1(k, 5) - 3)
        a, b = row_of[Y[j]], row_of[Y[i]]
        y = channel(a) + offset
        xs, ys = [X[j], X[j]], [bottom(a), y]
        for r in (a + 1):(b - 1)
            g = gutters[r][argmin(abs.(gutters[r] .- X[i]))]
            push!(xs, g, g); push!(ys, y, channel(r) + offset)
            y = channel(r) + offset
        end
        push!(xs, X[i], X[i]); push!(ys, y, Y[i] + 0.5)
        plot!(plt, xs, ys; color=colors[mod1(k, length(colors))], linewidth=1.2, label="")
    end
    return plt
end

# Each node is a unit box with a caption below it.
function _render_poset(poset, labels, drawings; per_row::Int, edges::Symbol, show_labels::Bool)
    X, Y, wrapped = _poset_layout(poset.covers; per_row)
    caption_gap = 0.45
    xlims = (minimum(X) - 0.8, maximum(X) + 0.8)
    ylims = (minimum(Y) - 0.6 - caption_gap, maximum(Y) + 0.7)
    plt = plot(; legend=:outerbottom, grid=false, axis=false, ticks=false, aspect_ratio=:equal, xlims, ylims, dpi=200,
               size=(max(500, round(Int, 175 * (xlims[2] - xlims[1]))), max(450, round(Int, 175 * (ylims[2] - ylims[1])))))

    for (bottom, top) in wrapped
        plot!(plt, Shape([xlims[1], xlims[2], xlims[2], xlims[1]],
                         [bottom - 0.5 - caption_gap, bottom - 0.5 - caption_gap, top + 0.6, top + 0.6]);
              fillcolor=:gray, fillalpha=0.08, linewidth=0, label="")
    end
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
