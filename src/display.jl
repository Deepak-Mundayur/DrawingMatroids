# REPL displays for DrawingMatroids result objects.
# Omit large fields such as systems, witness data,
# maximal-degeneration objects, complete histories, matrices, and plots.

_plural(n::Integer, singular::AbstractString, plural::AbstractString=singular * "s") =
    n == 1 ? singular : plural

function _matroid_description(M)
    try
        return "rank-$(Oscar.rank(M)) matroid on $(length(M)) elements"
    catch
        return string(typeof(M))
    end
end

function _parallel_class_count(reduction::DrawingReduction)
    return count(cls -> length(cls) > 1, reduction.parallel_classes)
end

function _displayed_point_count(reduction::DrawingReduction)
    return length(reduction.parallel_classes)
end

function _optional_diagnostic(diagnostics::AbstractDict, key)
    value = get(diagnostics, key, nothing)
    return value
end

function _format_seconds(value)
    value isa Real || return nothing
    value < 1e-3 && return @sprintf("%.3g ms", 1e3 * value)
    value < 60 && return @sprintf("%.3f s", value)
    minutes = floor(Int, value / 60)
    seconds = value - 60 * minutes
    return @sprintf("%d min %.1f s", minutes, seconds)
end

function _format_number(value::Real)
    value == 0 && return "0"
    magnitude = abs(value)
    if magnitude < 1e-3 || magnitude >= 1e4
        return @sprintf("%.3e", value)
    end
    return @sprintf("%.6g", value)
end

function Base.show(io::IO, initialization::InitializationResult)
    print(
        io,
        "InitializationResult(strategy=:", initialization.strategy,
        ", separating_bases=", length(initialization.separating_bases),
        ", maximal_degenerations=", length(initialization.maximal_degenerations),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", initialization::InitializationResult)
    diagnostics = initialization.diagnostics

    println(io, "InitializationResult")
    println(io, "  strategy:              :", initialization.strategy)
    println(io, "  state variables:       ", length(initialization.state))
    println(io, "  separating bases:      ", length(initialization.separating_bases))
    println(io, "  maximal degenerations: ", length(initialization.maximal_degenerations))

    component_index = initialization.component_index
    !isnothing(component_index) && println(io, "  component index:       ", component_index)

    attempt = _optional_diagnostic(diagnostics, :attempt)
    isnothing(attempt) && (attempt = _optional_diagnostic(diagnostics, :initializer_attempt))
    !isnothing(attempt) && println(io, "  successful attempt:    ", attempt)

    remi_seconds = _optional_diagnostic(diagnostics, :remi_seconds)
    formatted_time = _format_seconds(remi_seconds)
    !isnothing(formatted_time) && println(io, "  Rémi time:             ", formatted_time)

    rabinowitsch_count = _optional_diagnostic(diagnostics, :rabinowitsch_count)
    total_basis_count = _optional_diagnostic(diagnostics, :total_basis_count)
    if !isnothing(rabinowitsch_count) && !isnothing(total_basis_count)
        ratio = total_basis_count == 0 ? 0.0 : 100 * rabinowitsch_count / total_basis_count
        println(
            io,
            "  Rabinowitsch minors:   ", rabinowitsch_count,
            " / ", total_basis_count,
            " (", @sprintf("%.1f", ratio), "%)",
        )
    end
end

function Base.show(io::IO, validation::DrawingValidation)
    print(
        io,
        "DrawingValidation(valid=", validation.valid,
        ", max_nonbasis_residual=", _format_number(validation.maximum_nonbasis_residual),
        ", min_basis_determinant=", _format_number(validation.minimum_basis_absolute_determinant),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", validation::DrawingValidation)
    status = validation.valid ? "valid" : "invalid"
    println(io, "DrawingValidation (", status, ")")
    println(
        io,
        "  maximum nonbasis residual: ",
        _format_number(validation.maximum_nonbasis_residual),
    )
    println(
        io,
        "  minimum basis determinant: ",
        _format_number(validation.minimum_basis_absolute_determinant),
    )
    println(io, "  failed nonbases:            ", length(validation.failed_nonbases))
    println(io, "  failed bases:               ", length(validation.failed_bases))
    println(io, "  coordinate collisions:      ", length(validation.collisions))
    !isnothing(validation.component_member) &&
        println(io, "  component membership:       ", validation.component_member)
end

function Base.show(io::IO, reduction::DrawingReduction)
    print(
        io,
        "DrawingReduction(displayed_points=", _displayed_point_count(reduction),
        ", parallel_classes=", _parallel_class_count(reduction),
        ", loops=", length(reduction.loops),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", reduction::DrawingReduction)
    n_displayed = _displayed_point_count(reduction)
    n_parallel = _parallel_class_count(reduction)
    n_loops = length(reduction.loops)

    println(io, "DrawingReduction")
    println(io, "  target:               ", _matroid_description(reduction.matroid))
    println(io, "  displayed points:     ", n_displayed)
    println(io, "  nontrivial parallels: ", n_parallel)
    println(io, "  loops:                ", n_loops)

    if n_parallel > 0
        groups = [
            "{" * join(string.(labels), ",") * "}"
            for labels in reduction.grouped_labels
            if length(labels) > 1
        ]
        println(io, "  grouped labels:       ", join(groups, ", "))
    end

    n_loops > 0 && println(io, "  loop labels:          ", join(string.(reduction.loop_labels), ", "))
end

function Base.show(io::IO, result::DrawingResult)
    print(
        io,
        "DrawingResult(valid=", result.validation.valid,
        ", strategy=:", result.initialization.strategy,
        ", target=\"", _matroid_description(result.target_matroid), "\"",
        ", frames=", length(result.history),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", result::DrawingResult)
    reduction = result.reduction
    initialization = result.initialization
    validation = result.validation
    diagnostics = initialization.diagnostics

    n_elements = try
        length(result.target_matroid)
    catch
        length(reduction.matroid)
    end
    n_displayed = _displayed_point_count(reduction)
    n_parallel = _parallel_class_count(reduction)
    n_loops = length(reduction.loops)
    n_frames = length(result.history)

    println(io, "DrawingResult")
    println(io, "  target:            ", _matroid_description(result.target_matroid))
    println(io, "  status:            ", validation.valid ? "✓ valid realization" : "✗ invalid realization")
    println(io, "  initializer:       :", initialization.strategy)
    println(
        io,
        "  elements:          ", n_elements,
        " (", n_displayed, " displayed ", _plural(n_displayed, "point"),
        ", ", n_parallel, " parallel ", _plural(n_parallel, "class", "classes"),
        ", ", n_loops, " ", _plural(n_loops, "loop"), ")",
    )
    println(io, "  dynamics states:   ", n_frames)

    remi_seconds = _optional_diagnostic(diagnostics, :remi_seconds)
    formatted_time = _format_seconds(remi_seconds)
    !isnothing(formatted_time) && println(io, "  Rémi time:         ", formatted_time)

    n_degenerations = length(initialization.maximal_degenerations)
    n_separating = length(initialization.separating_bases)
    if n_degenerations > 0 || n_separating > 0
        println(io, "  degenerations:     ", n_degenerations)
        println(io, "  separating bases:  ", n_separating)
    end

    total_basis_count = _optional_diagnostic(diagnostics, :total_basis_count)
    if !isnothing(total_basis_count) && total_basis_count > 0
        println(
            io,
            "  Rabinowitsch use:  ", n_separating, " / ", total_basis_count,
            " basis minors (", @sprintf("%.1f", 100 * n_separating / total_basis_count), "%)",
        )
    end

    println(
        io,
        "  max NB residual:   ",
        _format_number(validation.maximum_nonbasis_residual),
    )
    println(
        io,
        "  min basis det:     ",
        _format_number(validation.minimum_basis_absolute_determinant),
    )

    if !isempty(validation.failed_nonbases) || !isempty(validation.failed_bases)
        println(io, "  failed nonbases:   ", length(validation.failed_nonbases))
        println(io, "  failed bases:      ", length(validation.failed_bases))
    end

    !isnothing(validation.component_member) &&
        println(io, "  component member:  ", validation.component_member)

    if !isnothing(result.image_filename) || !isnothing(result.animation_filename)
        println(io, "  output:")
        !isnothing(result.image_filename) && println(io, "    image:           ", result.image_filename)
        !isnothing(result.animation_filename) && println(io, "    animation:       ", result.animation_filename)
    end

    println(io, "  inspect: inspect(result) for more details")
end

function Base.show(io::IO, result::DrawingCollectionResult)
    print(io, "DrawingCollectionResult(drawings=", length(result.results), ")")
end

function Base.show(io::IO, ::MIME"text/plain", result::DrawingCollectionResult)
    println(io, "DrawingCollectionResult")
    println(io, "  drawings: ", length(result.results))
    valid_count = count(item -> item.validation.valid, result.results)
    println(io, "  valid:    ", valid_count, " / ", length(result.results))
    !isnothing(result.filename) && println(io, "  grid:     ", result.filename)

    for (index, item) in enumerate(result.results)
        println(
            io,
            "  [", index, "] ",
            _matroid_description(item.target_matroid),
            " — ",
            item.validation.valid ? "valid" : "invalid",
            " — :", item.initialization.strategy,
        )
    end
end

"""
    inspect(result::DrawingResult)

Print the user-relevant fields available on a drawing result.
The function does not print large arrays or systems themselves.
"""
function inspect(io::IO, result::DrawingResult)
    println(io, "DrawingResult inspection")
    println(io)
    println(io, "Primary result data")
    println(io, "  validation          result.validation")
    println(io, "  final coordinates   result.final_state")
    println(io, "  trajectory          result.history")
    println(io, "  plot                result.plot")
    println(io)
    println(io, "Initialization")
    println(io, "  summary             result.initialization")
    println(io, "  polynomial system   result.initialization.system")
    println(io, "  separating bases    result.initialization.separating_bases")
    println(io, "  degenerations       result.initialization.maximal_degenerations")
    println(io, "  diagnostics         result.initialization.diagnostics")
    println(io)
    println(io, "Realization")
    println(io, "  realization matrix  result.validation.matrix")
    println(io, "  target matroid      result.target_matroid")
    println(io, "  original source     result.source")
    println(io)
    println(io, "Combinatorial reduction")
    println(io, "  reduction           result.reduction")
    println(io, "  simple matroid      result.reduction.simple_matroid")
    println(io, "  parallel classes    result.reduction.parallel_classes")
    println(io, "  grouped labels      result.reduction.grouped_labels")
    println(io, "  loops               result.reduction.loop_labels")
    println(io)
    println(io, "Files")
    println(io, "  image               result.image_filename")
    println(io, "  animation           result.animation_filename")

    return nothing
end

inspect(result::DrawingResult) = inspect(stdout, result)