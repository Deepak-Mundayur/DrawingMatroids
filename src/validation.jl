# Map a drawing matrix `A` to a point in the variables of the parent space's system.
# Geometric coordinates come from `RealizationSpaces.coordinates` (ambient: vec(A);
# quotient: chart coordinates). Auxiliary (Rabinowitsch) variables enter the system
# affinely, F(x, u) = J_u * u + F(x, 0), so they are solved for directly.
function _component_coordinate_point(A::Matrix{Float64}, C; tol::Real)
    space = RealizationSpaces.parent(C)
    # A drawing of a lower-rank component is padded back to the space's matrix size.
    A = [A; zeros(RealizationSpaces.rank(space) - size(A, 1), size(A, 2))]
    F = RealizationSpaces.system(space)
    x = ComplexF64.(RealizationSpaces.coordinates(space, A; tol=tol, check=true))
    k = length(RealizationSpaces.auxiliary_variables(space))
    k == 0 && return x

    x0 = [x; zeros(ComplexF64, k)]
    J_u = HomotopyContinuation.jacobian(F, x0)[:, (end - k + 1):end]
    point = [x; J_u \ -F(x0)]
    all(isfinite, point) || return nothing
    maximum(abs, F(point)) <= sqrt(tol) * max(1.0, maximum(abs, point)) || return nothing
    return point
end

function _component_membership(A::Matrix{Float64}, C; tol::Real)
    isnothing(C) && return nothing
    point = _component_coordinate_point(A, C; tol=tol)
    isnothing(point) && return nothing
    return first(HomotopyContinuation.membership([point], RealizationSpaces.witness_set(C)))
end

# RealizationSpaces' realization diagnostics of the drawing, plus (when `component` is
# given) membership of the drawing in that component.
function compute_drawing_validation(state::AbstractVector, reduction; tol::Real=1e-7, component=nothing)
    A = coordinates_to_matrix(state[1:coordinate_count(reduction)], reduction)
    diagnostics = RealizationSpaces.matroid_realization_diagnostics(A, reduction.matroid; tol=tol)
    component_member = _component_membership(A, component; tol=tol)
    valid = diagnostics.valid && (isnothing(component) || component_member === true)
    return merge(diagnostics, (; valid, component_member))
end

# Generic matroid of C, checked against its witness points: each witness point is a generic
# point of C, so the nonbases are exactly the sets whose minors vanish at every witness point.
# RealizationSpaces computes the stored matroid from a single sample, which numerical error
# can occasionally turn into extra nonbases; that step is rerun for C until it agrees.
function _component_matroid!(C; tol::Real, attempts::Int=10)
    points = [RealizationSpaces.witness_point(C; point=i)
              for i in eachindex(HomotopyContinuation.solutions(RealizationSpaces.witness_set(C)))]
    r = maximum(A -> RealizationSpaces.numerical_rank(A; rtol=tol), points)
    witness_nonbases = Set(reduce(intersect, (RealizationSpaces.numerical_nonbases(A, r; rtol=tol) for A in points)))
    agrees(M) = !isnothing(M) && Oscar.rank(M) == r && Set(Oscar.nonbases(M)) == witness_nonbases

    for _ in 1:attempts
        agrees(RealizationSpaces.matroid(C)) && return RealizationSpaces.matroid(C)
        RealizationSpaces.compute_component_matroid!(C; tol=1e-9)
    end
    agrees(RealizationSpaces.matroid(C)) && return RealizationSpaces.matroid(C)
    error("RealizationSpaces' generic matroid of the component still disagrees with its witness points " *
          "after $attempts recomputations, so the component is not drawn.")
end
