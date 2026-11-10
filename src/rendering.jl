_group_label(labels) = length(labels) == 1 ? string(first(labels)) : "{" * join(string.(labels), ",") * "}"

# A label of two or more characters is set smaller so it fits inside its marker.
_marker_fontsize(label, fontsize) = length(label) == 1 ? fontsize : round(Int, 0.75 * fontsize)

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

# Offset of the point p from its nearest point on the segment ab.
function _segment_offset(px, py, ax, ay, bx, by)
    dx, dy = bx - ax, by - ay
    t = clamp(((px - ax) * dx + (py - ay) * dy) / (dx^2 + dy^2 + 1e-300), 0, 1)
    return px - ax - t * dx, py - ay - t * dy
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

    for flat in reduction.lines
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
                annotate!(plt, xs[i], ys[i], text(label, _marker_fontsize(label, fontsize), :white, :center))
            end
        end
    end

    loop_x, loop_y = _loop_positions(length(loops), bounds)
    if !isempty(loop_x)
        scatter!(plt, loop_x, loop_y; markershape=:rect, markersize=markersize + 1,
                 markercolor=:orange, markerstrokecolor=:orange, label=legend ? "loops" : "")
        if show_labels
            for i in eachindex(loop_x)
                annotate!(plt, loop_x[i], loop_y[i], text(string(loops[i]), _marker_fontsize(string(loops[i]), fontsize), :black, :center))
            end
        end
    end
    return plt
end

# A drawing as a plot: the box `bounds`, room above it for loops, no axes.
function render_drawing(state, reduction; bounds, framestyle, show_boundary, show_labels)
    (xmin, xmax), (ymin, ymax) = bounds
    top = isempty(Oscar.loops(reduction.matroid)) ? 0.5 : 1.5
    plt = plot(; xlim=(xmin - 0.5, xmax + 0.5), ylim=(ymin - 0.5, ymax + top), aspect_ratio=:equal,
               legend=:outerbottom, framestyle, xticks=false, yticks=false)
    return _draw_on!(plt, state[1:2length(reduction.parallel_classes)], reduction; bounds, show_boundary, show_labels)
end

# A gif of the optimizer's path, one frame per state.
function save_drawing_animation(history, reduction; filename, fps, style...)
    animation = Animation()
    foreach(state -> frame(animation, render_drawing(state, reduction; style...)), history)
    return gif(animation, filename; fps)
end
