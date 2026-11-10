# The result of draw_matroid; the REPL shows a short summary.
struct MatroidDrawing
    matrix::Matrix{Float64}
    plot::Plot
    seed::Int
    start::Symbol
end

Base.show(io::IO, d::MatroidDrawing) = print(io, "MatroidDrawing(", join(size(d.matrix), "×"), ", seed=", d.seed, ")")

function Base.show(io::IO, ::MIME"text/plain", d::MatroidDrawing)
    r, n = size(d.matrix)
    println(io, "Drawing of a rank-$r matroid on $n elements from start=$(repr(d.start)), seed=$(d.seed)")
    println(io, "  fields: matrix ($(r)×$n), plot, seed, start")
    println(io, "  matrix, rounded to 3 digits:")
    Base.print_matrix(IOContext(io, :compact => true), round.(d.matrix; digits=3), "    ")
end

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

Draw a matroid of rank at most three and return a `MatroidDrawing` with fields `matrix` (a
real realization matrix of `M`), `plot`, `seed` (it reproduces this drawing) and `start`
(the start that found it).

`start` chooses the initial point:

- `:maximal_degenerations`: maximal degenerations plus a greedy reduced Rabinowitsch
  transversal.
- `:circle`: points on a circle, pulled onto the drawing variety by the optimizer.
- `:exact`: a real point of Oscar's exact `realization_space(M; char=0)`, when no equations
  remain after Oscar's elimination or the space is finite.
- `:nid`: a real sample of a component of `RealizationSpaces.realization_space(M)`.
  The decomposition is the slow step; to draw several times, keep the component and
  call `draw_matroid(C)` on it instead.
- an r×n realization matrix of `M`.

If a start finds no drawing, `draw_matroid` continues with the next one in the order
`:maximal_degenerations`, `:circle`, `:exact`, `:nid`. The exact computation can also show
that `M` has no real realization. Before `:nid` it asks for confirmation, unless
`autorun=true`; outside an interactive session it then stops instead. A matrix start is
tried alone, and a component start (`draw_matroid(C)`) continues only with `:exact`.

Up to `attempts` attempts are made, all determined by `seed`. Without `seed` a random
one is chosen; `draw_matroid(M; seed=result.seed)` reproduces a drawing exactly, and
calling again without `seed` tries different random choices.

Keywords: `tol` (realization tolerance), `bounds`, `autorun`, `optimizer` (settings that
replace those in `DrawingMatroids.OPTIMIZER`, e.g. `optimizer=(k_line=10.0,)`), and the output
settings `filename`, `animate`, `animation_filename`, `fps`, `framestyle`, `show_boundary`,
`show_labels`.
"""
draw_matroid(M::Oscar.Matroid; kwargs...) = _draw_matroid(M, nothing; kwargs...)

# Optimizer settings: the step `dt`, `max_steps`, and the strengths of point repulsion `k_p`,
# the bounding box `k_wall`, the angle between lines `k_torque`, and of pushing points away from
# lines they are not on `k_line`.
const OPTIMIZER = (dt=0.01, max_steps=300, k_p=10.0, k_wall=50.0, k_torque=0.5, k_line=5.0)

# `component`: when given, the drawing must also lie on this RealizationComponent.
function _draw_matroid(
    M::Oscar.Matroid,
    component;
    start = :maximal_degenerations,
    attempts::Int = 5,
    seed::Union{Nothing,Integer} = nothing,
    tol::Real = 1e-7,
    bounds = ((-5.0, 5.0), (-5.0, 5.0)),
    autorun::Bool = false,
    optimizer::NamedTuple = (;),

    filename::Union{Nothing,AbstractString} = "matroid.png",
    animate::Bool = false,
    animation_filename::AbstractString = "matroid.gif",
    fps::Int = 15,
    framestyle::Symbol = :box,
    show_boundary::Bool = true,
    show_labels::Bool = true,
)
    attempts >= 1 || throw(ArgumentError("attempts must be positive."))
    keys(optimizer) ⊆ keys(OPTIMIZER) || throw(ArgumentError("optimizer settings are $(keys(OPTIMIZER))."))
    settings = merge(OPTIMIZER, optimizer)
    starts = (:maximal_degenerations, :circle, :exact, :nid)
    start in starts || start isa Union{AbstractMatrix, RealizationSpaces.RealizationComponent} || throw(ArgumentError(
        "Unsupported start=$(repr(start)). Use :maximal_degenerations, :circle, :exact, :nid, a RealizationComponent, " *
        "or a realization matrix."
    ))
    rank_value = Oscar.rank(M)
    rank_value <= 3 || error("DrawingMatroids only draws matroids of rank at most three (got rank $rank_value).")

    # Without nonbases the realization space is irreducible, so a generic circle works.
    start === :nid && isempty(Oscar.nonbases(M)) && (start = :circle)
    seed = Int(something(seed, rand(Random.RandomDevice(), UInt32)))
    style = (; bounds, framestyle, show_boundary, show_labels)

    chain = start isa Symbol ? starts[findfirst(==(start), starts):end] : start isa AbstractMatrix ? (start,) : (start, :exact)

    started = time()
    return _with_seed(seed) do
        reduction = drawing_reduction(M)
        cannot_draw(why) = error("DrawingMatroids cannot draw $(something(component, M)): $why.")
        last_error = nothing
        for s in chain
            s === :nid && s !== first(chain) && !autorun && !_confirm_nid() && break
            println("[DrawingMatroids] start=$(repr(_start_label(s)))")
            initial, member_of = s, component
            if s === :exact
                initial = exact_real_realization(M)
                initial isa AbstractString && cannot_draw(initial)
                if isnothing(initial)
                    last_error = ErrorException("start=:exact found no start: Oscar's realization space has equations left in positive dimension.")
                    println("[DrawingMatroids] start=:exact failed: Oscar's realization space has equations left in positive dimension")
                    continue
                end
            elseif s === :nid
                components = RealizationSpaces.components(RealizationSpaces.realization_space(M))
                isempty(components) && cannot_draw("it is not realizable over any field of characteristic 0")
                initial = member_of = first(components)
            end
            final_state, history = try
                _find_drawing(reduction, initial; attempts, seed, tol, bounds, settings, component=member_of)
            catch err
                err isa InterruptException && rethrow()
                last_error = err
                println("[DrawingMatroids] start=$(repr(_start_label(s))) failed: $(_first_line(err))")
                continue
            end
            println("[DrawingMatroids] Drew $(length(reduction.parallel_classes)) points from start=$(repr(_start_label(s))) in " *
                    "$(round(time() - started; digits=1)) s, seed=$seed" * (isnothing(filename) ? "" : ", saved $filename"))
            plt = render_drawing(final_state, reduction; style...)
            isnothing(filename) || savefig(plt, String(filename))
            animate && save_drawing_animation(history, reduction; filename=animation_filename, fps, style...)
            return MatroidDrawing(coordinates_to_matrix(final_state, reduction), plt, seed, _start_label(s))
        end
        throw(last_error)
    end
end

_start_label(s) = s isa Symbol ? s : s isa AbstractMatrix ? :matrix : :component

# The first line of an error message, at most 160 characters.
_first_line(err) = first(first(split(sprint(showerror, err), '\n')), 160)

# Ask before start=:nid, which decomposes the realization space and can take long.
function _confirm_nid()
    isinteractive() || (println("[DrawingMatroids] Skipping start=:nid; pass autorun=true to run it without asking."); return false)
    print("[DrawingMatroids] No drawing found yet. Try start=:nid (decomposes the realization space, can be slow)? [y/N] ")
    return lowercase(strip(readline())) in ("y", "yes")
end

# Try up to `attempts` times to reach an exact drawing; return (final_state, history).
# A matrix start is validated once; :circle, :maximal_degenerations and component starts get a
# fresh point on each attempt (a rotated, jittered circle, for :maximal_degenerations projected
# onto the reduced Rabinowitsch system first).
function _find_drawing(reduction, start; attempts, seed, tol, bounds, settings, component)
    simple = reduction.simple_matroid
    point_count = length(reduction.parallel_classes)

    # Any distinct points on a line realize a simple matroid of rank at most two.
    if Oscar.rank(simple) < 3
        coordinates = start isa AbstractMatrix ? realization_matrix_to_coordinates(start, reduction; bounds, real_tol=tol) :
                      low_rank_default_coordinates(point_count, bounds)
        validation = compute_drawing_validation(coordinates, reduction; tol, component)
        validation.valid || error("The drawing does not realize the matroid: $validation")
        return coordinates, [coordinates]
    end

    if start isa AbstractMatrix
        fixed = realization_matrix_to_coordinates(start, reduction; bounds, real_tol=tol)
        validation = compute_drawing_validation(fixed, reduction; tol)
        validation.valid || throw(ArgumentError("The start point does not realize the matroid: $validation"))
    end
    system = drawing_system(simple)
    start === :maximal_degenerations && (rabinowitsch = rabinowitsch_system(simple))
    # At most 5% of the box width per optimizer step (the step is dt times the force).
    max_force = 0.05 * (bounds[1][2] - bounds[1][1]) / settings.dt

    last_validation = last_error = nothing
    for attempt in 1:attempts
        try
            p0, guess = if start isa AbstractMatrix
                fixed, nothing
            elseif start in (:circle, :maximal_degenerations)
                guess = circle_guess(point_count, bounds; phase=attempt == 1 ? 0.0 : 2pi * rand(), jitter=0.015 * (attempt - 1))
                start === :circle ? (nothing, guess) : (project_rabinowitsch(rabinowitsch, guess, reduction; tol), nothing)
            else
                sampled = start_coordinates(start, reduction; bounds, real_tol=tol)
                last_validation = compute_drawing_validation(sampled, reduction; tol)
                last_validation.valid || continue
                sampled, nothing
            end

            # Each further attempt takes smaller steps under a larger force cap.
            vector_field = drawing_force_field(point_count, bounds, reduction.lines; settings.k_p, settings.k_wall,
                settings.k_torque, settings.k_line, max_force=max_force * 2.0^(attempt - 1))
            final_state, history = constrained_optimize(system, vector_field; p0, guess, dt=settings.dt / 2.0^(attempt - 1),
                settings.max_steps, corrector=moore_penrose_corrector, tol=1e-10, max_iters=100)

            last_validation = compute_drawing_validation(final_state, reduction; tol, component)
            if !last_validation.valid
                # The optimizer can slide into a degeneration (points merging); keep the last
                # state along its path that still realizes the matroid.
                k = findlast(state -> compute_drawing_validation(state, reduction; tol).valid, history)
                if !isnothing(k)
                    final_state, history = history[k], history[1:k]
                    last_validation = compute_drawing_validation(final_state, reduction; tol, component)
                end
            end
            if last_validation.valid
                history = [Float64.(state) for state in history]
                spread = _readable_frame(_ipopt_polish(Float64.(final_state), reduction, bounds; tol), reduction, bounds; tol)
                return spread, push!(history, spread)
            end
        catch err
            err isa InterruptException && rethrow()
            last_error = err
        end
    end

    error("No drawing found after $attempts attempts (seed=$seed); try another seed. Last error: " *
          "$(isnothing(last_error) ? "none" : _first_line(last_error)). Last validation: $last_validation." *
          (start isa RealizationSpaces.RealizationComponent ? " The target is the matroid stored on the component, " *
           "which RealizationSpaces may have computed wrongly." : ""))
end

"""
    draw_matroid(C::RealizationComponent; check_component_membership=true, ...)

Draw the generic matroid of the irreducible component `C` (for example a component of
`RealizationSpaces.realization_space(M)` or `nonbasis_variety(M)`), starting from a real
sample of `C` and, by default, checking that the final drawing lies on `C`. Returns a
`MatroidDrawing`; other keywords are those of `draw_matroid(M)`.
"""
draw_matroid(C::RealizationSpaces.RealizationComponent; kwargs...) = _draw_component(C, _component_matroid!(C; tol=1e-7); kwargs...)

# draw_matroid(C) with the generic matroid M of C already checked by _component_matroid!.
function _draw_component(C, M; check_component_membership::Bool=true, kwargs...)
    haskey(kwargs, :start) && throw(ArgumentError("draw_matroid(C) always starts from a sample of C."))
    return _draw_matroid(M, check_component_membership ? C : nothing; kwargs..., start=C)
end
