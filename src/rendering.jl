_group_label(labels) = length(labels) == 1 ? string(first(labels)) : "{" * join(string.(labels), ",") * "}"

# Segment between the two outermost points of a line (flat), ordered along the line.
function _line_segment(coords::AbstractVector, flat::Vector{Int})
    xs = [coords[2i - 1] for i in flat]
    ys = [coords[2i] for i in flat]
    far = argmax(hypot.(xs .- xs[1], ys .- ys[1]))   # gives the direction of the line
    dx, dy = xs[far] - xs[1], ys[far] - ys[1]
    t =(xs .- xs[1]) .* dx .+ (ys .- ys[1]) .* dy
    ends = [argmin(t), argmax(t)]
    return xs[ends], ys[ends]
end

# Axis-aligned boxes (xmin, xmax, ymin, ymax): distance to a point, overlap area, overflow of `bounds`.
_box_distance(b, (px, py)) = hypot(max(b[1] - px, 0, px - b[2]), max(b[3] - py, 0, py - b[4]))
_box_overlap(b, q) = max(0, min(b[2], q[2]) - max(b[1], q[1])) * max(0, min(b[4], q[4]) - max(b[3], q[3]))
_box_outside(b, ((xmin, xmax), (ymin, ymax))) = max(xmin - b[1], 0) + max(b[2] - xmax, 0) + max(ymin - b[3], 0) + max(b[4] - ymax, 0)

# Where to put the label of a parallel class at (x, y): right, left, above or below the point.
# Prefer a spot inside `bounds`, then one not overlapping labels already placed, then the one
# farthest from the other points. Returns (x, y, halign, valign, box).
function _label_spot(x, y, w, h, offset, others, placed, bounds)
    (xmin, xmax), (ymin, ymax) = bounds
    spots = [(x + offset, y, :left, :vcenter, (x + offset, x + offset + w, y - h / 2, y + h / 2)),
             (x - offset, y, :right, :vcenter, (x - offset - w, x - offset, y - h / 2, y + h / 2)),
             (x, y + offset, :center, :bottom, (x - w / 2, x + w / 2, y + offset, y + offset + h)),
             (x, y - offset, :center, :top, (x - w / 2, x + w / 2, y - offset - h, y - offset))]
    score(s) = (-_box_outside(s[5], bounds), -sum((_box_overlap(s[5], q) for q in placed); init=0.0),
                minimum((_box_distance(s[5], p) for p in others); init=Inf))
    return argmax(score, spots)
end

# Label placement for every parallel class with more than one element: class index => spot.
function _label_layout(coords, reduction, bounds)
    width = bounds[1][2] - bounds[1][1]
    offset, char_width = 0.03 * width, 0.025 * width
    k = length(reduction.parallel_classes)
    placed, layout = NTuple{4,Float64}[], Dict{Int,Any}()
    for i in 1:k
        labels = reduction.parallel_classes[i]
        length(labels) == 1 && continue
        label = _group_label(labels)
        others = [(coords[2j - 1], coords[2j]) for j in 1:k if j != i]
        spot = _label_spot(coords[2i - 1], coords[2i], length(label) * char_width, 1.6 * char_width,
                           offset, others, placed, bounds)
        push!(placed, spot[5])
        layout[i] = spot
    end
    return layout
end

# Readability thresholds, as fractions of the box width.
const _MIN_POINT_GAP = 0.10        # between two points
const _MIN_LINE_GAP = 0.06         # between a point and a drawn line it is not on

# Offset of the point p from its nearest point on the segment ab.
function _segment_offset(px, py, ax, ay, bx, by)
    dx, dy = bx - ax, by - ay
    t = clamp(((px - ax) * dx + (py - ay) * dy) / (dx^2 + dy^2 + 1e-300), 0, 1)
    return px - ax - t * dx, py - ay - t * dy
end

# Readability of a drawing inside `bounds`: the smallest of (closest two points) / _MIN_POINT_GAP
# and (closest point to a drawn line it is not on) / _MIN_LINE_GAP, in box widths; 0 if a label
# leaves the box or overlaps a point or another label. A drawing is readable when this is >= 1.
function _readability(coords, reduction, bounds)
    width = bounds[1][2] - bounds[1][1]
    k = length(reduction.parallel_classes)
    point(i) = (coords[2i - 1], coords[2i])
    score = Inf
    for i in 1:k, j in (i + 1):k
        score = min(score, hypot((point(i) .- point(j))...) / (_MIN_POINT_GAP * width))
    end
    for flat in rank_two_flats(reduction.simple_matroid)
        (ax, bx), (ay, by) = _line_segment(coords, flat)
        for p in 1:k
            p in flat || (score = min(score, hypot(_segment_offset(point(p)..., ax, ay, bx, by)...) / (_MIN_LINE_GAP * width)))
        end
    end
    boxes = [spot[5] for spot in values(_label_layout(coords, reduction, bounds))]
    marker = 0.015 * width
    for (a, b) in enumerate(boxes)
        _box_outside(b, bounds) > 0 && return 0.0
        any(_box_overlap(b, q) > 0 for q in boxes[(a + 1):end]) && return 0.0
    end
    # A label's own point sits 0.03 widths from its box; no point may come within a marker radius.
    any(_box_distance(b, point(i)) < marker for b in boxes for i in 1:k) && return 0.0
    return score
end

# A projective transformation keeps every line and every non-collinearity, so it turns a
# drawing into another drawing of the same matroid; it acts by an invertible matrix on the left,
# under which every component of a realization space is invariant, so it also stays on the
# same component. A readable drawing is kept; otherwise the frame is chosen to maximize
# `_readability` (random search, then local refinement), and used if it still validates.
function _readable_frame(coords, reduction, bounds; tol, candidates=200, refinements=400)
    k = length(reduction.parallel_classes)
    k >= 2 || return coords
    H = [ones(1, k); reshape(coords[1:2k], 2, k)]
    function framed(G)
        B = G * H
        all(abs.(B[1, :]) .> 1e-6) || return nothing
        return _fit_coordinates_to_bounds(vec(B[2:3, :] ./ B[1:1, :]), bounds)
    end
    score(c) = isnothing(c) ? -Inf : _readability(c, reduction, bounds)
    original = score(framed(Matrix(1.0I, 3, 3)))
    original >= 1 && return coords
    best_G, best = Matrix(1.0I, 3, 3), original
    for _ in 1:candidates
        G = randn(3, 3)
        s = score(framed(G))
        s > best && ((best_G, best) = (G, s))
    end
    step = 0.3
    for _ in 1:refinements
        best >= 1.5 && break
        G = best_G .+ step * sqrt(sum(abs2, best_G) / 9) .* randn(3, 3)
        s = score(framed(G))
        s > best ? ((best_G, best) = (G, s)) : (step *= 0.995)
    end
    best > original || return coords
    c = framed(best_G)
    return compute_drawing_validation(c, reduction; tol).valid ? c : coords
end


# Loops are shown in a row just above the drawing box.
function _loop_positions(loop_count::Int, bounds)
    loop_count == 0 && return Float64[], Float64[]
    (xmin, xmax), (ymin, ymax) = bounds
    xs = loop_count == 1 ? [(xmin + xmax) / 2] : collect(range(xmin, xmax; length=loop_count))
    ys = fill(ymax + 0.075 * (ymax - ymin), loop_count)
    return xs, ys
end

# Draw the points, lines, labels and loops of a drawing into `plt`, inside `bounds`.
function _draw_on!(plt, coords, reduction; bounds, show_boundary::Bool, show_labels::Bool,
                   markersize=8, fontsize=9, legend::Bool=true)
    point_count = length(reduction.parallel_classes)
    (xmin, xmax), (ymin, ymax) = bounds
    loops = sort!(Int.(collect(Oscar.loops(reduction.matroid))))

    if show_boundary
        plot!(plt, [xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin];
              linewidth=1.5, linestyle=:dash, color=:red, label="")
    end

    for flat in rank_two_flats(reduction.simple_matroid)
        line_x, line_y = _line_segment(coords, flat)
        plot!(plt, line_x, line_y; linewidth=1.5, alpha=0.55, color=:gray, label="")
    end

    xs = [coords[2i - 1] for i in 1:point_count]
    ys = [coords[2i] for i in 1:point_count]
    scatter!(plt, xs, ys; markershape=:circle, markersize=markersize,
             markercolor=:blue, markerstrokecolor=:blue,
             label=legend ? "points / parallel classes" : "")

    if show_labels
        layout = _label_layout(coords, reduction, bounds)
        for i in 1:point_count
            label = _group_label(reduction.parallel_classes[i])
            if haskey(layout, i)
                lx, ly, halign, valign, _ = layout[i]
                annotate!(plt, lx, ly, text(label, fontsize, :black, halign, valign))
            else
                annotate!(plt, xs[i], ys[i], text(label, fontsize, :white, :center))
            end
        end
    end

    loop_x, loop_y = _loop_positions(length(loops), bounds)
    if !isempty(loop_x)
        scatter!(plt, loop_x, loop_y; markershape=:rect, markersize=markersize + 1,
                 markercolor=:orange, markerstrokecolor=:orange, label=legend ? "loops" : "")
        if show_labels
            for i in eachindex(loop_x)
                annotate!(plt, loop_x[i], loop_y[i], text(string(loops[i]), fontsize, :black, :center))
            end
        end
    end
    return plt
end

function render_drawing(
    state::AbstractVector,
    reduction;
    bounds,
    framestyle,
    show_boundary::Bool,
    show_labels::Bool,
)
    (xmin, xmax), (ymin, ymax) = bounds
    extra_top = isempty(Oscar.loops(reduction.matroid)) ? 0.5 : 1.5
    plt = plot(
        xlim=(xmin - 0.5, xmax + 0.5),
        ylim=(ymin - 0.5, ymax + extra_top),
        aspect_ratio=:equal,
        legend=:outerbottom,
        title="",
        framestyle=framestyle,
        xticks = false,
        yticks= false
    )
    coords = state[1:(2 * length(reduction.parallel_classes))]
    return _draw_on!(plt, coords, reduction; bounds, show_boundary, show_labels)
end

function save_static_drawing(
    state,
    reduction;
    filename::Union{Nothing,AbstractString},
    style...,
)
    plt = render_drawing(state, reduction; style...)
    if !isnothing(filename)
        savefig(plt, String(filename))
    end
    return plt
end

function save_drawing_animation(
    history::AbstractVector,
    reduction;
    filename::AbstractString,
    fps::Int,
    style...,
)
    animation = Animation()
    for state in history
        frame(animation, render_drawing(state, reduction; style...))
    end

    gif(animation, filename; fps=fps)
    return String(filename)
end