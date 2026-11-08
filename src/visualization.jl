# Run `f` with the global RNG seeded by `seed`, then restore the caller's RNG state.
# RealizationSpaces samples with the global RNG, so this covers its randomness too.
function _with_seed(f, seed::Integer)
    saved = copy(Random.default_rng())
    Random.seed!(seed)
    try
        return f()
    finally
        copy!(Random.default_rng(), saved)
    end
end

"""
    draw_matroid(M::Matroid; start=:maximal_degenerations, seed=nothing, ...)

Draw a rank-two or rank-three matroid and return `(matrix, plot, seed)`: a real
realization matrix of `M`, the plot, and the seed that reproduces this drawing.

`start` chooses the initial point:

- `:maximal_degenerations`: maximal degenerations plus a greedy reduced Rabinowitsch
  transversal.
- `:circle`: points on a circle, pulled onto the drawing variety by the optimizer.
- `:nid`: a real sample of a component of `RealizationSpaces.realization_space(M)`.
  The decomposition is the slow step; to draw several times, keep the component and
  call `draw_component` on it instead.
- an r×n realization matrix of `M`.

Up to `attempts` attempts are made, all determined by `seed`. Without `seed` a random
one is chosen; `draw_matroid(M; seed=result.seed)` reproduces a drawing exactly, and
calling again without `seed` tries different random choices.

Keywords: `tol` (realization tolerance), `bounds`, the optimizer settings `dt`,
`max_steps`, `k_p` (point repulsion), `k_wall` (bounding box), `k_torque` (angle between
lines), and the output settings `filename`, `animate`, `animation_filename`, `fps`,
`framestyle`, `show_boundary`, `show_labels`.
"""
function draw_matroid(
    M::Oscar.Matroid;
    start = :maximal_degenerations,
    attempts::Int = 5,
    seed::Union{Nothing,Integer} = nothing,
    tol::Real = 1e-7,
    bounds = ((-5.0, 5.0), (-5.0, 5.0)),

    dt::Real = 0.01,
    max_steps::Int = 300,
    k_p::Real = 10.0,
    k_wall::Real = 50.0,
    k_torque::Real = 0.5,

    filename::Union{Nothing,AbstractString} = "matroid.png",
    animate::Bool = false,
    animation_filename::AbstractString = "matroid.gif",
    fps::Int = 15,
    framestyle::Symbol = :box,
    show_boundary::Bool = true,
    show_labels::Bool = true,

    # Set by draw_component.
    component = nothing,
    check_component_membership::Bool = false,
)
    attempts >= 1 || throw(ArgumentError("attempts must be positive."))
    start in (:maximal_degenerations, :circle, :nid) || start isa Union{AbstractMatrix, RealizationSpaces.RealizationComponent} || throw(ArgumentError(
        "Unsupported start=$(repr(start)). Use :maximal_degenerations, :circle, :nid, a RealizationComponent, " *
        "or a realization matrix."
    ))
    rank_value = Oscar.rank(M)
    rank_value in 2:3 || error("DrawingMatroids only draws matroids of rank two or three (got rank $rank_value).")
    Oscar.is_realizable(M; char=0) || error(
        "DrawingMatroids cannot draw $(something(component, M)): it is not realizable over any field of characteristic 0."
    )

    # Without nonbases the realization space is irreducible, so a generic circle works.
    start === :nid && isempty(Oscar.nonbases(M)) && (start = :circle)
    seed = Int(something(seed, rand(Random.RandomDevice(), UInt32)))
    style = (; bounds, framestyle, show_boundary, show_labels)

    return _with_seed(seed) do
        initial, member_of = start, check_component_membership ? component : nothing
        if start === :nid
            initial = member_of = first(RealizationSpaces.components(RealizationSpaces.realization_space(M)))
        end

        reduction = drawing_reduction(M)
        final_state, history = _find_drawing(reduction, initial;
            attempts, seed, tol, bounds, dt, max_steps, k_p, k_wall, k_torque, component=member_of)

        plt = save_static_drawing(final_state, reduction; filename, style...)
        animate && save_drawing_animation(history, reduction; filename=animation_filename, fps, style...)
        (matrix=coordinates_to_matrix(final_state, reduction), plot=plt, seed=seed)
    end
end

# Try up to `attempts` times to reach an exact drawing; return (final_state, history).
# Fixed starts (maximal degenerations, user input) are computed once and validated;
# :circle and component starts get a fresh random point on each attempt.
function _find_drawing(reduction, start;
    attempts, seed, tol, bounds, dt, max_steps, k_p, k_wall, k_torque, component)

    simple = reduction.simple_matroid
    point_count = length(reduction.parallel_classes)
    rank_three = Oscar.rank(simple) == 3

    fixed_start = if !rank_three && start isa Symbol
        rank_two_default_coordinates(point_count, bounds)
    elseif start === :maximal_degenerations
        maximal_degeneration_start(reduction; bounds, attempts, tol)
    elseif start isa AbstractMatrix
        realization_matrix_to_coordinates(start, reduction; bounds, real_tol=tol)
    end

    if !isnothing(fixed_start)
        validation = compute_drawing_validation(fixed_start, reduction; tol)
        validation.valid || throw(ArgumentError("The start point does not realize the matroid: $validation"))
    end

    if rank_three
        system = drawing_system(simple)
        vector_field = drawing_force_field(point_count, bounds, rank_two_flats(simple); k_p, k_wall, k_torque)
    end

    last_validation = last_error = nothing
    for attempt in 1:attempts
        try
            p0, guess = if !isnothing(fixed_start)
                fixed_start, nothing
            elseif start === :circle
                nothing, circle_guess(point_count, bounds;
                    phase=attempt == 1 ? 0.0 : 2pi * rand(), jitter=0.015 * (attempt - 1))
            else
                sampled = start_coordinates(start, reduction; bounds, real_tol=tol)
                last_validation = compute_drawing_validation(sampled, reduction; tol)
                last_validation.valid || continue
                sampled, nothing
            end

            final_state, history = !rank_three ? (p0, [p0]) : constrained_optimize(system, vector_field;
                p0, guess, dt=dt / 2.0^(attempt - 1), max_steps,
                corrector=moore_penrose_corrector, tol=1e-10, max_iters=100)

            last_validation = compute_drawing_validation(final_state, reduction; tol, component)
            last_validation.valid && return Float64.(final_state), [Float64.(state) for state in history]
        catch err
            err isa InterruptException && rethrow()
            last_error = err
        end
    end

    error(
        "No exact drawing was found after $attempts attempts with seed=$seed; " *
        "rerun without `seed` (or with another seed) to try different random choices. " *
        "Last validation: $last_validation. Last error: $(isnothing(last_error) ? "none" : sprint(showerror, last_error))." *
        (start isa RealizationSpaces.RealizationComponent ?
            " The drawing targets the matroid with nonbases $(Oscar.nonbases(reduction.matroid)), the generic " *
            "matroid stored on the component (RealizationSpaces.matroid(C)). If that is not the component's " *
            "matroid, RealizationSpaces may have stored a wrong one due to numerical errors." : "")
    )
end

"""
    draw_component(C::RealizationComponent; start=C, check_component_membership=true, ...)

Draw the generic matroid of `C`, starting from a real sample of `C` and, by default,
checking that the final drawing lies on `C`. Returns `(matrix, plot, seed)`; other
keywords are those of `draw_matroid`.
"""
function draw_component(
    C::RealizationSpaces.RealizationComponent;
    start=C,
    check_component_membership::Bool=true,
    kwargs...,
)
    M = _component_matroid!(C; tol=1e-7)
    return draw_matroid(M;
        start=start === :nid ? C : start,
        component=C, check_component_membership, kwargs...,
    )
end
