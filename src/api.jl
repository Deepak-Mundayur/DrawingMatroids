"""
    run_and_animate_collinear_system(n_pts, collinear_sets; kwargs...)

Compatibility wrapper around Deepak's original circle-initialized workflow.
It always builds the Constrained_Optimization system with
`hard_non_collinearity_relns=false`.
"""
function run_and_animate_collinear_system(
    n_pts::Int,
    collinear_sets;
    bounds=((-5.0, 5.0), (-5.0, 5.0)),
    dt::Real=0.01,
    max_steps::Int=300,
    k_p::Real=10.0,
    k_wall::Real=50.0,
    guess=nothing,
    filename::AbstractString="example.gif",
    fps::Int=15,
    framestyle=:none,
)
    triples = NTuple{3,Int}[
        (Int(t[1]), Int(t[2]), Int(t[3]))
        for t in collinear_sets
    ]
    system = matroid_collinearity_system(n_pts, triples)
    initial_guess = isnothing(guess) ? circle_guess(n_pts, bounds) : Float64.(collect(guess))
    field = Constrained_Optimization.make_bounded_repelling_force(
        n_pts,
        bounds;
        k_p=k_p,
        k_wall=k_wall,
    )

    final_point, history = constrained_optimize(
        system,
        field;
        guess=initial_guess,
        dt=dt,
        max_steps=max_steps,
        corrector=Constrained_Optimization.moore_penrose_corrector,
    )

    animation = Animation()
    (xmin, xmax), (ymin, ymax) = bounds

    for (step, state) in enumerate(history)
        xs = [state[2 * i - 1] for i in 1:n_pts]
        ys = [state[2 * i] for i in 1:n_pts]
        plt = plot(
            xlim=(xmin - 1.0, xmax + 1.0),
            ylim=(ymin - 1.0, ymax + 1.0),
            aspect_ratio=:equal,
            legend=false,
            title="DrawingMatroids.jl (step $step)",
            framestyle=framestyle,
        )

        for triple in triples
            endpoints = _line_endpoints(state, collect(triple), bounds)
            isnothing(endpoints) && continue
            line_x, line_y = endpoints
            plot!(plt, line_x, line_y; alpha=0.4, linewidth=1.5)
        end

        scatter!(plt, xs, ys; markersize=8)
        frame(animation, plt)
    end

    gif(animation, filename; fps=fps)
    return final_point, history
end

function run_and_animate_collinear_system(n_pts::Int, bounds, collinear_sets; kwargs...)
    return run_and_animate_collinear_system(
        n_pts,
        collinear_sets;
        bounds=bounds,
        kwargs...,
    )
end
