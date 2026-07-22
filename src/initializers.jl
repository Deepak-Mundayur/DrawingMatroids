function _target_component_matroid(C; tol::Real=1e-7)
    if isnothing(RealizationSpaces.matroid(C))
        RealizationSpaces.compute_component_matroid!(C; tol=tol)
    end
    return RealizationSpaces.matroid(C)
end

function _matching_components(RS, target; tol::Real=1e-7)
    cps = RealizationSpaces.components(RS)
    matches = Tuple{Int,Any}[]

    for (index, C) in enumerate(cps)
        if isnothing(RealizationSpaces.matroid(C))
            RealizationSpaces.compute_component_matroid!(C; tol=tol)
        end

        CM = RealizationSpaces.matroid(C)

        RealizationSpaces.same_labeled_matroid(CM, target) &&
            push!(matches, (index, C))
    end

    return matches
end

function nid_initial_coordinates(target, reduction::DrawingReduction; bounds=((-5.0, 5.0), (-5.0, 5.0)), real_tol::Real=1e-7, matroid_tol::Real=1e-7, rng::Random.AbstractRNG=Random.default_rng())
    if isempty(collect(Oscar.nonbases(target)))
        coordinates = Oscar.rank(target) == 2 ? 
        rank_two_default_coordinates(length(reduction.parallel_classes), bounds) : circle_guess(length(reduction.parallel_classes), bounds)
        return (coordinates=coordinates, component=nothing, component_index=nothing, realization_space=nothing, diagnostics=Dict{Symbol,Any}(:ambient_component => true))
    end

    Oscar.is_realizable(target; char=0) || error("The matroid is not realizable over a field of characteristic zero, and therefore is not realizable over C.")

    RS = RealizationSpaces.nonbasis_variety(target; compute_NID=true, compute_orbits=false)
    matches = _matching_components(RS, target; tol=matroid_tol)

    isempty(matches) && error("Oscar.is_realizable(target; char=0) reports that the matroid is realizable in characteristic zero, but the numerical irreducible decomposition did not contain a component whose computed generic labeled matroid equals the target. This indicates an NID, witness, tolerance, or component-matroid reconstruction failure.")

    last_error = nothing

    for (index, C) in matches
        try
            A = RealizationSpaces.sample(C; real_point=true)
            Ar = Matrix{Float64}(real.(A))
            RealizationSpaces.satisfies_matroid(Ar, target; tol=matroid_tol) || continue

            coordinates = realization_matrix_to_coordinates(Ar, reduction; bounds=bounds, real_tol=real_tol, rng=rng)

            return (coordinates=coordinates, component=C, component_index=index, realization_space=RS, diagnostics=Dict{Symbol,Any}(:matching_component_count => length(matches), :real_source => :realization_spaces_sample))
        catch err
            last_error = err
        end
    end

    error("The matroid is realizable in characteristic zero and has at least one matching complex realization component, but no valid numerically real sample suitable for drawing was obtained. Last error: $(repr(last_error)).")
end

function component_initial_coordinates(C, reduction::DrawingReduction; bounds=((-5.0, 5.0), (-5.0, 5.0)), real_tol::Real=1e-7, matroid_tol::Real=1e-7, rng::Random.AbstractRNG=Random.default_rng())
    target = reduction.relabeled.matroid

    A = RealizationSpaces.sample(C; real_point=true)
    Ar = Matrix{Float64}(real.(A))

    RealizationSpaces.satisfies_matroid(Ar, target; tol=matroid_tol) || error(
        "The real point sampled from the requested RealizationComponent does not realize the expected labeled target matroid."
    )

    return realization_matrix_to_coordinates(Ar, reduction; bounds=bounds, real_tol=real_tol, rng=rng)
end

function _rabinowitsch_initial_point(
    reduction::DrawingReduction,
    separating_bases::Vector{Vector{Int}};
    bounds,
    phase::Real=0.0,
    jitter::Real=0.0,
    attempts::Int=5,
    projection_tol::Real=1e-10,
    projection_max_iters::Int=100,
    matroid_tol::Real=1e-7,
    rng::Random.AbstractRNG=Random.default_rng(),
)
    attempts >= 1 || throw(ArgumentError("initializer attempts must be positive."))

    simple = reduction.simple_matroid
    point_count = length(simple)
    coordinate_length = 2 * point_count
    initializer_system = drawing_system(simple; separating_bases=separating_bases)

    last_error = nothing
    last_validation = nothing

    for attempt in 1:attempts
        current_phase = attempt == 1 ? phase : 2pi * rand(rng)
        current_jitter = attempt == 1 ? jitter : max(jitter, 0.015 * (attempt - 1))

        try
            affine_guess = circle_guess(point_count, bounds; phase=current_phase, jitter=current_jitter, rng=rng,
            )
            augmented_guess = append_rabinowitsch_coordinates(affine_guess, separating_bases,
            )

            println(
                "[DrawingMatroids] Projecting reduced Rabinowitsch initializer " *
                "(attempt $attempt/$attempts, $(length(separating_bases)) separating bases)...",
            )

            augmented_point, success = Constrained_Optimization.project_onto_manifold(
                initializer_system, augmented_guess;
                tol=projection_tol, max_iters=projection_max_iters,
            )

            success || begin
                last_error = ErrorException(
                    "Reduced Rabinowitsch projection did not converge on attempt $attempt."
                )
                continue
            end

            coordinates = Float64.(augmented_point[1:coordinate_length])
            validation = validate_drawing(
                coordinates,
                reduction;
                tol=matroid_tol,
                collision_tol=sqrt(Float64(matroid_tol)),
            )
            last_validation = validation

            if validation.valid
                return (
                    coordinates=coordinates,
                    initializer_system=initializer_system,
                    attempt=attempt,
                    validation=validation,
                )
            end

            last_error = ErrorException(
                "Projected initializer did not realize the target matroid: " *
                string(validation_summary(validation)),
            )
        catch err
            last_error = err
        end
    end

    validation_text = isnothing(last_validation) ?
        "no initializer validation was available" :
        string(validation_summary(last_validation))

    error(
        "Could not obtain an exact real initializer from the Rabinowitsch system " *
        "after $attempts attempts. Last validation: $validation_text. " *
        "Last error: $(repr(last_error))."
    )
end

function prepare_rank_three_initialization(
    reduction::DrawingReduction;
    strategy::Symbol=:maximal_degenerations,
    start_point=nothing,
    source_component=nothing,
    bounds=((-5.0, 5.0), (-5.0, 5.0)),
    phase::Real=0.0,
    jitter::Real=0.0,
    initializer_attempts::Int=5,
    projection_tol::Real=1e-10,
    projection_max_iters::Int=100,
    remi_preprocess::Bool=true,
    remi_verbosity::Int=0,
    remi_threads::Int=Base.Threads.nthreads(),
    real_tol::Real=1e-7,
    matroid_tol::Real=1e-7,
    rng::Random.AbstractRNG=Random.default_rng(),
)
    simple = reduction.simple_matroid
    target = reduction.relabeled.matroid
    n_points = length(simple)

    separating_bases = Vector{Vector{Int}}()
    degenerations = Any[]
    component_index = nothing
    diagnostics = Dict{Symbol,Any}()
    p0 = nothing
    guess = nothing

    system = drawing_system(simple)
    initializer_system = system

    if strategy in (:maximal_degenerations, :remi)
        println("[DrawingMatroids] Computing Rémi maximal degenerations with $(min(remi_threads, Base.Threads.nthreads())) Julia thread(s)...")

        remi_start = time()
        separating_bases, degenerations = maximal_degeneration_basis_transversal(
            simple;
            preprocess=remi_preprocess, verbosity=remi_verbosity, n_procs=remi_threads
        )
        diagnostics[:remi_seconds] = time() - remi_start

        println("[DrawingMatroids] Rémi finished: $(length(degenerations)) maximal degenerations and $(length(separating_bases)) separating bases.")

        initializer_data = _rabinowitsch_initial_point(reduction, separating_bases;
            bounds=bounds, phase=phase, jitter=jitter, attempts=initializer_attempts,
            projection_tol=projection_tol, projection_max_iters=projection_max_iters,
            matroid_tol=matroid_tol, rng=rng
        )
        p0 = initializer_data.coordinates
        initializer_system = initializer_data.initializer_system

        total_basis_count = length(collect(Oscar.bases(simple)))
        diagnostics[:initializer_attempt] = initializer_data.attempt
        diagnostics[:initializer_source] = :remi_reduced_rabinowitsch
        diagnostics[:rabinowitsch_count] = length(separating_bases)
        diagnostics[:total_basis_count] = total_basis_count
        diagnostics[:rabinowitsch_reduction_ratio] =
            total_basis_count == 0 ? 0.0 : length(separating_bases) / total_basis_count
        diagnostics[:maximal_degeneration_count] = length(degenerations)
        diagnostics[:dynamics_uses_rabinowitsch] = false

    elseif strategy == :all_bases
        separating_bases = sort(
            [sort(Int.(collect(B))) for B in Oscar.bases(simple)];
            by=B -> Tuple(B)
        )

        initializer_data = _rabinowitsch_initial_point(reduction, separating_bases; bounds=bounds, phase=phase, jitter=jitter,
            attempts=initializer_attempts,
            projection_tol=projection_tol, projection_max_iters=projection_max_iters,
            matroid_tol=matroid_tol, rng=rng
        )
        p0 = initializer_data.coordinates
        initializer_system = initializer_data.initializer_system

        diagnostics[:initializer_attempt] = initializer_data.attempt
        diagnostics[:initializer_source] = :all_bases_rabinowitsch
        diagnostics[:rabinowitsch_count] = length(separating_bases)
        diagnostics[:total_basis_count] = length(separating_bases)
        diagnostics[:rabinowitsch_reduction_ratio] = 1.0
        diagnostics[:baseline] = :all_bases
        diagnostics[:dynamics_uses_rabinowitsch] = false

    elseif strategy == :deepak
        guess = circle_guess(n_points, bounds; phase=phase, jitter=jitter, rng=rng)
        diagnostics[:hard_non_collinearity_relns] = false
        diagnostics[:initializer_source] = :deepak_circle

    elseif strategy == :nid
        data = nid_initial_coordinates(target, reduction; bounds=bounds, real_tol=real_tol, matroid_tol=matroid_tol, rng=rng
        )
        p0 = data.coordinates
        component_index = data.component_index
        merge!(diagnostics, data.diagnostics)
        diagnostics[:realization_space] = data.realization_space
        diagnostics[:component] = data.component
        diagnostics[:initializer_source] = :nid

    elseif strategy == :user
        p0 = user_input_to_coordinates(start_point, reduction; bounds=bounds, real_tol=real_tol, rng=rng
        )
        diagnostics[:initializer_source] = :user

    elseif strategy == :component
        isnothing(source_component) && throw(ArgumentError(
            "start_point_strategy=:component requires a RealizationComponent source."
        ))
        p0 = component_initial_coordinates(source_component, reduction; bounds=bounds, real_tol=real_tol, matroid_tol=matroid_tol, rng=rng
        )
        diagnostics[:component] = source_component
        diagnostics[:initializer_source] = :component

    else
        throw(ArgumentError("Unknown start_point_strategy=$strategy. Use :maximal_degenerations, :all_bases, :nid, :deepak, :user, or :component."
        ))
    end

    return (system=system, initializer_system=initializer_system, p0=p0,
        guess=guess, separating_bases=separating_bases, maximal_degenerations=degenerations,
        component_index=component_index, diagnostics=diagnostics,
    )
end
