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

function _loop_positions(loop_count::Int, bounds)
    loop_count == 0 && return Float64[], Float64[]
    (xmin, xmax), (_, ymax) = bounds
    xs = loop_count == 1 ? [(xmin + xmax) / 2] : collect(range(xmin, xmax; length=loop_count))
    ys = fill(ymax + 0.75, loop_count)
    return xs, ys
end

function render_drawing(
    state::AbstractVector,
    reduction;
    bounds,
    framestyle,
    show_boundary::Bool,
    show_labels::Bool,
)
    simple = reduction.simple_matroid
    point_count = length(reduction.parallel_classes)
    coords = state[1:(2 * point_count)]
    (xmin, xmax), (ymin, ymax) = bounds
    loops = sort!(Int.(collect(Oscar.loops(reduction.matroid))))

    extra_top = isempty(loops) ? 0.5 : 1.5

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

    if show_boundary
        plot!(plt, [xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin];
              linewidth=1.5, linestyle=:dash, color=:red, label="")
    end

    for flat in rank_two_flats(simple)
        line_x, line_y = _line_segment(coords, flat)
        plot!(plt, line_x, line_y; linewidth=1.5, alpha=0.55, color=:gray, label="")
    end

    xs = [coords[2i - 1] for i in 1:point_count]
    ys = [coords[2i] for i in 1:point_count]

    scatter!(plt, xs, ys; markershape=:circle, markersize=8,
             markercolor=:blue, markerstrokecolor=:blue,
             label="points / parallel classes")

    if show_labels
        parallel_label_offset = 0.03 * (xmax - xmin)
        for i in 1:point_count
            labels = reduction.parallel_classes[i]
            label = _group_label(labels)
            if length(labels) == 1
                annotate!(plt, xs[i], ys[i], text(label, 9, :white, :center))
            else
                annotate!(plt, xs[i] + parallel_label_offset, ys[i], text(label, 9, :black, :left))
            end
        end
    end

    loop_x, loop_y = _loop_positions(length(loops), bounds)

    if !isempty(loop_x)
        scatter!(plt, loop_x, loop_y; markershape=:rect, markersize=9,
                 markercolor=:orange, markerstrokecolor=:orange,
                 label="loops")

        if show_labels
            for i in eachindex(loop_x)
                annotate!(plt, loop_x[i], loop_y[i],
                          text(string(loops[i]), 9, :black, :center))
            end
        end
    end

    return plt
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