function _coordinate_collisions(coords::AbstractVector, point_count::Int; tol::Real=1e-6)
    collisions = Tuple{Int,Int}[]
    for i in 1:point_count, j in (i + 1):point_count
        dx = coords[2 * i - 1] - coords[2 * j - 1]
        dy = coords[2 * i] - coords[2 * j]
        hypot(dx, dy) <= tol && push!(collisions, (i, j))
    end
    return collisions
end

function _component_coordinate_point(A::Matrix{Float64}, C; tol::Real=1e-10)
    parent_space = RealizationSpaces.parent(C)
    parent_system = RealizationSpaces.system(parent_space)
    variable_count = length(variables(parent_system))

    if RealizationSpaces.is_quotient_space(parent_space)
        coordinates = ComplexF64.(RealizationSpaces.quotient_coordinates(
            parent_space, A; tol=tol, check=true))
        chart = RealizationSpaces.quotient_chart(parent_space)

        if chart.has_rabinowitsch_variable
            open_conditions = RealizationSpaces._nontrivial_inequations(chart.inequations)
            product_value = prod(
                RealizationSpaces._evaluate_oscar_polynomial(f, coordinates)
                for f in open_conditions
            )
            abs(product_value) > tol || return nothing
            push!(coordinates, inv(product_value))
        end

        length(coordinates) == variable_count || return nothing
        return coordinates
    end

    variable_count == length(A) || return nothing
    return vec(ComplexF64.(A))
end

function _component_membership(A::Matrix{Float64}, C)
    isnothing(C) && return nothing

    try
        point = _component_coordinate_point(A, C)
        isnothing(point) && return nothing

        return first(HomotopyContinuation.membership(
            [point],
            RealizationSpaces.witness_set(C),
        ))
    catch
        return nothing
    end
end

function validate_drawing(state::AbstractVector, reduction::DrawingReduction;
    tol::Real=1e-7,
    collision_tol::Real=1e-6,
    component=nothing,
)
    coordinate_length = coordinate_count(reduction)
    length(state) >= coordinate_length || throw(ArgumentError(
        "Drawing state has length $(length(state)), but at least $coordinate_length coordinates are required."
    ))

    coords = state[1:coordinate_length]
    A = coordinates_to_matrix(coords, reduction)

    minors = RealizationSpaces.matroid_realization_diagnostics(A, reduction.matroid; tol=tol)

    collisions = _coordinate_collisions(coords, length(reduction.parallel_classes); tol=collision_tol)
    component_member = _component_membership(A, component)

    valid = minors.valid
    if !isnothing(component)
        valid &= (component_member === true)
    end

    return DrawingValidation(valid, A,
        minors.max_nonbasis_residual,
        minors.min_basis_absolute_determinant,
        minors.failed_nonbases,
        minors.failed_bases,
        collisions,
        component_member,
    )
end

function validation_summary(validation::DrawingValidation)
    return (
        valid=validation.valid,
        maximum_nonbasis_residual=validation.maximum_nonbasis_residual,
        minimum_basis_absolute_determinant=validation.minimum_basis_absolute_determinant,
        failed_nonbases=validation.failed_nonbases,
        failed_bases=validation.failed_bases,
        collisions=validation.collisions,
        component_member=validation.component_member,
    )
end