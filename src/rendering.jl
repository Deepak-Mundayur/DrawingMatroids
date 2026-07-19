_group_label(labels) = length(labels) == 1 ? string(first(labels)) : "{" * join(string.(labels), ",") * "}"

function _line_endpoints(coords::AbstractVector, flat::Vector{Int}, bounds; tol::Real=1e-10)
    (xmin, xmax), (ymin, ymax) = bounds
    length(flat) >= 2 || return nothing

    p = flat[1]
    q = nothing
    xp, yp = coords[2p - 1], coords[2p]

    for candidate in flat[2:end]
        xc, yc = coords[2 * candidate - 1], coords[2 * candidate]
        if hypot(xc - xp, yc - yp) > tol
            q = candidate
            break
        end
    end

    isnothing(q) && return nothing
    xq, yq = coords[2q - 1], coords[2q]
    dx, dy = xq - xp, yq - yp

    if abs(dx) > tol
        x_values = [xmin - 1.0, xmax + 1.0]
        slope = dy / dx
        intercept = yp - slope * xp
        y_values = slope .* x_values .+ intercept
        return x_values, y_values
    end

    return [xp, xp], [ymin - 1.0, ymax + 1.0]
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
    reduction::DrawingReduction;
    bounds=((-5.0, 5.0), (-5.0, 5.0)),
    title::AbstractString="DrawingMatroids.jl",
    framestyle=:box,
    show_boundary::Bool=true,
    show_labels::Bool=true,
)
    simple = reduction.simple_matroid
    point_count = length(reduction.parallel_classes)
    coords = state[1:(2 * point_count)]
    (xmin, xmax), (ymin, ymax) = bounds

    extra_top = isempty(reduction.loops) ? 0.5 : 1.5
    plt = plot(
        xlim=(xmin - 0.5, xmax + 0.5),
        ylim=(ymin - 0.5, ymax + extra_top),
        aspect_ratio=:equal,
        legend=:outertopright,
        title=title,
        framestyle=framestyle,
    )

    if show_boundary
        plot!(
            plt,
            [xmin, xmax, xmax, xmin, xmin],
            [ymin, ymin, ymax, ymax, ymin];
            linewidth=1.5,
            linestyle=:dash,
            color=:red,
            label="",
        )
    end

    flats = rank_two_flats(simple)
    line_label_used = false
    for flat in flats
        endpoints = _line_endpoints(coords, flat, bounds)
        isnothing(endpoints) && continue
        line_x, line_y = endpoints
        plot!(
            plt,
            line_x,
            line_y;
            linewidth=1.5,
            alpha=0.55,
            color=:gray,
            label=line_label_used ? "" : "rank-2 flats",
        )
        line_label_used = true
    end

    xs = [coords[2 * i - 1] for i in 1:point_count]
    ys = [coords[2 * i] for i in 1:point_count]
    scatter!(
        plt,
        xs,
        ys;
        markershape=:circle,
        markersize=8,
        markercolor=:blue,
        markerstrokecolor=:blue,
        label="points / parallel classes",
    )

    if show_labels
        # Singleton labels are centered inside their blue marker. Parallel-class
        # labels are placed just to the right so grouped labels do not cover the
        # marker itself.
        parallel_label_offset = 0.03 * (xmax - xmin)

        for i in 1:point_count
            labels = reduction.grouped_labels[i]
            label = _group_label(labels)

            if length(labels) == 1
                annotate!(
                    plt,
                    xs[i],
                    ys[i],
                    text(label, 9, :white, :center),
                )
            else
                annotate!(
                    plt,
                    xs[i] + parallel_label_offset,
                    ys[i],
                    text(label, 9, :black, :left),
                )
            end
        end
    end

    loop_x, loop_y = _loop_positions(length(reduction.loops), bounds)
    if !isempty(loop_x)
        scatter!(
            plt,
            loop_x,
            loop_y;
            markershape=:rect,
            markersize=9,
            markercolor=:orange,
            markerstrokecolor=:orange,
            label="", # Loops are intentionally omitted from the legend.
        )

        if show_labels
            for i in eachindex(loop_x)
                annotate!(plt, loop_x[i], loop_y[i], text(string(reduction.loop_labels[i]), 9, :left))
            end
        end
    end

    return plt
end

function save_static_drawing(
    state,
    reduction;
    filename::Union{Nothing,AbstractString}=nothing,
    kwargs...,
)
    plt = render_drawing(state, reduction; kwargs...)
    if !isnothing(filename)
        savefig(plt, String(filename))
    end
    return plt
end

function save_drawing_animation(
    history::AbstractVector,
    reduction::DrawingReduction;
    filename::AbstractString="drawing.gif",
    fps::Int=15,
    bounds=((-5.0, 5.0), (-5.0, 5.0)),
    title::AbstractString="DrawingMatroids.jl",
    framestyle=:box,
    show_boundary::Bool=true,
    show_labels::Bool=true,
)
    animation = Animation()

    for (step, state) in enumerate(history)
        plt = render_drawing(
            state,
            reduction;
            bounds=bounds,
            title="$title (step $step)",
            framestyle=framestyle,
            show_boundary=show_boundary,
            show_labels=show_labels,
        )
        frame(animation, plt)
    end

    gif(animation, filename; fps=fps)
    return String(filename)
end