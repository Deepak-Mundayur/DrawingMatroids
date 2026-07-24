function _rank_two_system(point_count::Int)
    vars = Variable[]
    for i in 1:point_count
        push!(vars, Variable("x$i"), Variable("y$i"))
    end
    return System([0.0 * vars[1]]; variables=vars)
end

function _validation_component(source_component, diagnostics)
    !isnothing(source_component) && return source_component
    return get(diagnostics, :component, nothing)
end

function _rank_two_initialization(
    reduction::DrawingReduction;
    strategy::Symbol, start_point,
    source_component, bounds,
    real_tol, matroid_tol, rng)

    simple = reduction.simple_matroid
    point_count = length(simple)
    diagnostics = Dict{Symbol,Any}()
    component_index = nothing

    coordinates = if strategy == :user
        user_input_to_coordinates(start_point, reduction; bounds=bounds, real_tol=real_tol, rng=rng)
    elseif strategy == :nid
        data = nid_initial_coordinates(reduction.matroid, reduction;
            bounds=bounds, real_tol=real_tol, matroid_tol=matroid_tol,
            rng=rng)
        component_index = data.component_index
        merge!(diagnostics, data.diagnostics)
        diagnostics[:component] = data.component
        data.coordinates
    elseif strategy == :component
        isnothing(source_component) && throw(ArgumentError(
            "start_point_strategy=:component requires a RealizationComponent."
        ))
        component_initial_coordinates(source_component, reduction;
            bounds=bounds, real_tol=real_tol, matroid_tol=matroid_tol, rng=rng)
    elseif strategy in (:maximal_degenerations, :remi, :deepak)
        diagnostics[:rank_two_closed_form] = true
        rank_two_default_coordinates(point_count, bounds)
    else
        throw(ArgumentError(
            "Unknown rank-two start_point_strategy=$strategy."
        ))
    end

    system = _rank_two_system(point_count)
    return InitializationResult(strategy, Float64.(coordinates),
        system, Vector{Vector{Int}}(), Any[], component_index, diagnostics)
end

function _drawing_result(source, target_matroid, reduction, initialization,
    final_state, history, validation;
    bounds, title, framestyle, show_boundary, show_labels, filename,
    animate, animation_filename, fps,
)
    plt = save_static_drawing(final_state, reduction;
        filename=filename, bounds=bounds, title=title,
        framestyle=framestyle, show_boundary=show_boundary, show_labels=show_labels)

    gif_file = nothing
    if animate
        gif_file = isnothing(animation_filename) ? "drawing.gif" : String(animation_filename)
        save_drawing_animation(history, reduction;
            filename=gif_file, fps=fps, bounds=bounds, title=title,
            framestyle=framestyle, show_boundary=show_boundary, show_labels=show_labels)
    end

    return DrawingResult(source, target_matroid, reduction, initialization,
        Float64.(final_state), [Float64.(state) for state in history],
        validation, plt, isnothing(filename) ? nothing : String(filename),
        gif_file,
    )
end

"""
    visualization(M::Matroid; start_point_strategy=:maximal_degenerations, ...)

Draw a rank-two or rank-three matroid. Rank zero and rank one will throw
an error. The supported start strategies are:

- `:maximal_degenerations` / `:remi`: Rémi's maximal degenerations plus a greedy
  reduced Rabinowitsch transversal.
- `:nid`: RealizationSpaces numerical irreducible decomposition.
- `:deepak`: Deepak's circle guess and Newton pull-in without Rémi or NID.
- `:user`: user-provided realization matrix or affine coordinates.
"""
function visualization(
    M::Oscar.Matroid;
    start_point_strategy::Symbol = :maximal_degenerations,
    start_point = nothing,

    # Used internally by the RealizationComponent and RealizationSpace methods.
    source_component = nothing,
    source_object = nothing,
    check_component_membership::Bool = false,

    bounds = ((-5.0, 5.0), (-5.0, 5.0)),
    attempts::Int = 5,
    seed = nothing,

    filename::Union{Nothing,AbstractString} = "matroid.png",
    animate::Bool = false,
    animation_filename::Union{Nothing,AbstractString} = "matroid.gif",

    options::DrawingOptions = DrawingOptions(),
)
    attempts >= 1 || throw(ArgumentError("attempts must be positive."))
    rng = isnothing(seed) ? Random.default_rng() : Random.MersenneTwister(seed)

    reduction = drawing_reduction(M)
    target = reduction.matroid
    rank_value = Oscar.rank(target)
    result_source = isnothing(source_object) ? M : source_object

    rank_value == 0 && error("DrawingMatroids does not draw rank-zero matroids.")
    rank_value == 1 && error("DrawingMatroids does not draw rank-one matroids.")
    rank_value <= 3 || error("DrawingMatroids supports matroids of rank at most three.")

    if rank_value == 2
        initialization = _rank_two_initialization(reduction;
            strategy=start_point_strategy, start_point=start_point,
            source_component=source_component, bounds=bounds,
            real_tol=options.real_tol, matroid_tol=options.validation_tol, rng=rng)
        final_state = initialization.state
        history = [copy(final_state)]
        component_for_validation = check_component_membership ?
            _validation_component(source_component, initialization.diagnostics) : nothing
        validation = validate_drawing(final_state, reduction;
            tol=options.validation_tol, collision_tol=options.collision_tol,
            component=component_for_validation)
        validation.valid || error(
            "The rank-two initializer did not produce the requested matroid: $(validation_summary(validation))."
        )

        return _drawing_result(
            result_source, target, reduction, initialization,
            final_state, history, validation;
            bounds=bounds, title=options.title, framestyle=options.framestyle,
            show_boundary=options.show_boundary, show_labels=options.show_labels,
            filename=filename, animate=animate,
            animation_filename=animation_filename, fps=options.fps,
        )
    end

    prepared = prepare_rank_three_initialization(reduction;
        strategy=start_point_strategy, start_point=start_point,
        source_component=source_component, bounds=bounds,
        phase=0.0, jitter=0.0,
        remi_preprocess=options.remi_preprocess, remi_verbosity=options.remi_verbosity,
        remi_threads=options.remi_threads, initializer_attempts=attempts,
        projection_tol=options.projection_tol, projection_max_iters=options.projection_max_iters,
        real_tol=options.real_tol, matroid_tol=options.validation_tol, rng=rng,
    )

    if prepared.p0 !== nothing
        initial_validation = validate_drawing(prepared.p0, reduction;
            tol=options.validation_tol, collision_tol=options.collision_tol,
        )
        initial_validation.valid || throw(ArgumentError(
            "The supplied initializer does not realize the requested matroid: $(validation_summary(initial_validation))."
        ))
    end

    point_count = length(reduction.parallel_classes)
    vector_field = make_bounded_repelling_force(
        point_count, bounds;
        k_p=options.k_p, k_wall=options.k_wall,
    )

    last_validation = nothing
    last_error = nothing

    for attempt in 1:attempts
        local p0 = prepared.p0
        local guess = prepared.guess
        local current_dt = options.dt / 2.0^(attempt - 1)

        # Rémi/all-bases/NID/user/component strategies already provide p0.
        # Only Deepak's comparison strategy needs a fresh projected guess.
        if start_point_strategy == :deepak && attempt > 1
            p0 = nothing
            guess = circle_guess(point_count, bounds; phase=2pi * rand(rng), jitter=0.015 * (attempt - 1), rng=rng)
        end

        try
            final_state, history = constrained_optimize(
                prepared.system, vector_field;
                p0=p0, guess=guess, dt=current_dt, max_steps=options.max_steps,
                corrector=moore_penrose_corrector, tol=options.projection_tol, max_iters=options.projection_max_iters,
            )

            component_for_validation = check_component_membership ?
                _validation_component(source_component, prepared.diagnostics) : nothing
            validation = validate_drawing(final_state, reduction;
                tol=options.validation_tol, collision_tol=options.collision_tol,
                component=component_for_validation,
            )
            last_validation = validation

            if validation.valid
                diagnostics = copy(prepared.diagnostics)
                diagnostics[:attempt] = attempt
                diagnostics[:dt_used] = current_dt
                diagnostics[:hard_non_collinearity_relns] = false

                initialization = InitializationResult(
                    start_point_strategy,
                    Float64.(first(history)),
                    prepared.initializer_system,
                    prepared.separating_bases,
                    prepared.maximal_degenerations,
                    prepared.component_index,
                    diagnostics,
                )

                return _drawing_result(result_source, target, reduction,
                    initialization, final_state, history, validation;
                    bounds=bounds, title=options.title, framestyle=options.framestyle,
                    show_boundary=options.show_boundary, show_labels=options.show_labels,
                    filename=filename, animate=animate,
                    animation_filename=animation_filename, fps=options.fps,
                )
            end
        catch err
            last_error = err
        end
    end

    diagnostic_text = isnothing(last_validation) ? "no final validation was available" :
        string(validation_summary(last_validation))
    error(
        "No exact drawing was found after $attempts attempts using strategy $start_point_strategy. " *
        "Last validation: $diagnostic_text. Last error: $(repr(last_error))."
    )
end

function visualization(
    C::RealizationSpaces.RealizationComponent;
    start_point_strategy::Symbol=:component,
    check_component_membership::Bool=true,
    kwargs...,
)
    target = _target_component_matroid(C)
    effective_strategy = start_point_strategy == :nid ? :component : start_point_strategy
    return visualization(target;
        start_point_strategy=effective_strategy,
        source_component=C, source_object=C,
        check_component_membership=check_component_membership,
        kwargs...,
    )
end

function _ensure_components!(RS)
    if isnothing(RealizationSpaces.stored_components(RS))
        RealizationSpaces.get_NID!(RS)
    end
    return RealizationSpaces.components(RS)
end

function _component_filename(filename, index::Int)
    isnothing(filename) && return nothing
    path = String(filename)
    stem, extension = splitext(path)
    return "$(stem)_component_$(index)$(extension)"
end

function _collection_grid(results::Vector{DrawingResult}; filename=nothing)
    isempty(results) && error("No component drawings were produced.")
    count = length(results)
    columns = ceil(Int, sqrt(count))
    rows = ceil(Int, count / columns)
    combined = plot([result.plot for result in results]...; layout=(rows, columns))
    !isnothing(filename) && savefig(combined, String(filename))
    return DrawingCollectionResult(results, combined, isnothing(filename) ? nothing : String(filename))
end

function visualization(
    RS::RealizationSpaces.RealizationSpace;
    all::Bool=false,
    component::Int=0,
    start_point_strategy::Symbol=:component,
    filename::Union{Nothing,AbstractString}=nothing,
    grid_filename::Union{Nothing,AbstractString}=nothing,
    kwargs...,
)
    target_if_matroid = RealizationSpaces.structure_kind(RS) == :matroid ?
        RealizationSpaces.matroid(RS) : nothing
    ambient_matroid_space = !isnothing(target_if_matroid) &&
        isempty(collect(Oscar.nonbases(target_if_matroid)))

    if all && ambient_matroid_space
        result = visualization(
            target_if_matroid;
            start_point_strategy=start_point_strategy == :component ? :deepak : start_point_strategy,
            filename=_component_filename(filename, 1),
            kwargs...,
        )
        return _collection_grid([result]; filename=grid_filename)
    end

    if component > 0 && ambient_matroid_space
        component == 1 || throw(BoundsError([1], component))
        return visualization(
            target_if_matroid;
            start_point_strategy=start_point_strategy == :component ? :deepak : start_point_strategy,
            filename=filename,
            kwargs...,
        )
    end

    if all
        cps = _ensure_components!(RS)
        results = DrawingResult[]
        for i in eachindex(cps)
            push!(results, visualization(
                cps[i];
                start_point_strategy=start_point_strategy,
                filename=_component_filename(filename, i),
                kwargs...,
            ))
        end
        return _collection_grid(results; filename=grid_filename)
    end

    if component > 0
        cps = _ensure_components!(RS)
        component in eachindex(cps) || throw(BoundsError(cps, component))
        return visualization(cps[component];
            start_point_strategy=start_point_strategy,
            filename=filename,
            kwargs...,
        )
    end

    if RealizationSpaces.structure_kind(RS) == :incidence
        error(
            "An incidence-defined RealizationSpace requires `component=i` or `all=true` so that a component matroid can be inferred."
        )
    end

    target = RealizationSpaces.matroid(RS)
    isnothing(target) && error("The RealizationSpace does not store a parent matroid.")

    if isempty(collect(Oscar.nonbases(target)))
        return visualization(target;
            start_point_strategy=start_point_strategy == :component ? :deepak : start_point_strategy,
            filename=filename, kwargs...,
        )
    end

    Oscar.is_realizable(target; char=0) || error(
        "The parent matroid of the RealizationSpace is not realizable over characteristic zero."
    )
    _ensure_components!(RS)
    matches = _matching_components(RS, target)
    isempty(matches) && error(
        "The numerical irreducible decomposition did not contain a component whose computed generic labeled matroid equals the characteristic-zero-realizable parent matroid."
    )
    selected_index, selected_component = first(matches)

    result = visualization(selected_component;
        start_point_strategy=start_point_strategy,
        filename=filename,
        kwargs...,
    )
    result.initialization.diagnostics[:selected_component_index] = selected_index
    return result
end

matroid_visualize(args...; kwargs...) = visualization(args...; kwargs...)

# function matroid_visualize(M)

#     return nothing
# end



"""
    run_and_animate_collinear_system(n_pts, bounds, collinear_sets; kwargs...)

Creates the matroid geometry, defines the force field, finds an equilibrium, and outputs a GIF.
"""
function run_and_animate_collinear_system(n_pts::Int, collinear_sets; 
                                          bounds = ((-5.0,5.0),(-5.0,5.0)),
                                          dt=0.01, max_steps=300, 
                                          k_p=10.0, k_wall=50.0,
                                          guess = nothing,
                                          filename="example.gif", fps=15,
                                          show_guess = false,
                                          framestyle = :none,
                                          hard_non_collinearity_relns = false,
                                          )
    
    (x_min, x_max), (y_min, y_max) = bounds

    
    F_sys = _matroid_collinearity_system(n_pts, collinear_sets; hard_non_collinearity_relns=hard_non_collinearity_relns)
    n_vars = length(variables(F_sys))
    
    V = make_bounded_repelling_force(n_pts, bounds; k_p=k_p, k_wall=k_wall)
    
    if guess == nothing
        guess = Float64[]
        
        radius = min(x_max - x_min, y_max - y_min) * 0.35
        center_x = (x_max + x_min) / 2.0
        center_y = (y_max + y_min) / 2.0
        
        for i in 1:n_pts
            # Distribute points evenly around a circle
            angle = 2 * pi * i / n_pts
            push!(guess, center_x + radius * cos(angle))
            push!(guess, center_y + radius * sin(angle))
        end
        
        # Fill remaining slack variables if using the hard-boundary matroid system
        num_slacks = n_vars - length(guess)
        if num_slacks > 0
            append!(guess, fill(0.5, num_slacks))
        end
    end    

    guess_xs = [guess[2*i - 1] for i in 1:n_pts]
    guess_ys = [guess[2*i] for i in 1:n_pts]
    
    final_point, history = constrained_optimize(F_sys, V; guess=guess, dt=dt, max_steps=max_steps, corrector = moore_penrose_corrector)
    anim = Animation()
    
    for (step_idx, state) in enumerate(history)
        xs = [state[2*i - 1] for i in 1:n_pts]
        ys = [state[2*i] for i in 1:n_pts]
        
        plt = plot(
            xlim=(x_min - 1.0, x_max + 1.0), 
            ylim=(y_min - 1.0, y_max + 1.0),
            aspect_ratio=:equal,
            legend=false,
            title="DrawingMatroids.jl",
            framestyle = framestyle
        )
        
        plot!(plt, [x_min, x_max, x_max, x_min, x_min], 
                   [y_min, y_min, y_max, y_max, y_min], 
                   color=:red, linewidth=2, linestyle=:dash)
        
        for (p1, p2, p3) in collinear_sets
            dx = xs[p2] - xs[p1]
            dy = ys[p2] - ys[p1]
            
            # Safety check: If p1 and p2 are completely crushed together 
            # (like in step 1), use p3 to calculate the line instead
            if abs(dx) < 1e-5 && abs(dy) < 1e-5
                dx = xs[p3] - xs[p1]
                dy = ys[p3] - ys[p1]
            end
            
            # If all 3 points are crushed into the exact same pixel, skip the line for this frame
            if abs(dx) < 1e-5 && abs(dy) < 1e-5
                continue
            end

            if abs(dx) > 1e-7
                slope = dy / dx
                intercept = ys[p1] - slope * xs[p1]
                line_xs = [x_min - 2.0, x_max + 2.0]
                line_ys = slope .* line_xs .+ intercept
                plot!(plt, line_xs, line_ys, color=:gray, alpha=0.4, linewidth=1.5)
            else
                # Handle perfectly vertical lines
                plot!(plt, [xs[p1], xs[p1]], [y_min - 2.0, y_max + 2.0], color=:gray, alpha=0.4, linewidth=1.5)
            end
        end

        # Scatter the initial guess in the background
        if show_guess
            scatter!(plt, guess_xs, guess_ys, color=:lightgray, markersize=8, alpha=0.5)
        end 

        # Scatter the active beads
        scatter!(plt, xs, ys, color=:blue, markersize=8)
        
        frame(anim, plt)
    end
    
    gif(anim, filename, fps=fps)
    println("Animation saved to: $filename")
    
    return final_point, history
end


function run_and_animate_collinear_system(n_pts::Int, bounds, collinear_sets; kwargs...)
    return run_and_animate_collinear_system(n_pts, collinear_sets; bounds = bounds, kwargs...)
end
