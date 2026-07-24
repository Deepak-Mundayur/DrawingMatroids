function collinearity_triples(M)
    Oscar.rank(M) == 3 || return NTuple{3,Int}[]

    triples = NTuple{3,Int}[]
    for NB in Oscar.nonbases(M)
        indices = sort(Int.(collect(NB)))
        length(indices) == 3 || continue
        push!(triples, (indices[1], indices[2], indices[3]))
    end

    sort!(triples)
    return unique(triples)
end

function _affine_minor_expression(vars_xy, basis::AbstractVector{<:Integer})
    length(basis) == 3 || throw(ArgumentError("Rank-three affine minors require three indices."))
    i, j, k = Int.(basis)

    xi, yi = vars_xy[2 * i - 1], vars_xy[2 * i]
    xj, yj = vars_xy[2 * j - 1], vars_xy[2 * j]
    xk, yk = vars_xy[2 * k - 1], vars_xy[2 * k]

    return (xj - xi) * (yk - yi) - (yj - yi) * (xk - xi)
end

function affine_minor_value(coords::AbstractVector, basis::AbstractVector{<:Integer})
    length(basis) == 3 || throw(ArgumentError("Rank-three affine minors require three indices."))
    i, j, k = Int.(basis)

    xi, yi = coords[2 * i - 1], coords[2 * i]
    xj, yj = coords[2 * j - 1], coords[2 * j]
    xk, yk = coords[2 * k - 1], coords[2 * k]

    return (xj - xi) * (yk - yi) - (yj - yi) * (xk - xi)
end

"""
Extract a HomotopyContinuation `System` from supported
Constrained_Optimization return formats.

Some Constrained_Optimization versions return a `System` directly, while
others return `(name, System)`.
"""
_as_hc_system(system::System) = system

function _as_hc_system(result::Tuple)
    found = System[]

    for item in result
        item isa System && push!(found, item)
    end

    length(found) == 1 || throw(ArgumentError(
        "Expected exactly one HomotopyContinuation System in the " *
        "Constrained_Optimization result, but found $(length(found))."
    ))

    return only(found)
end

function _as_hc_system(result)
    throw(ArgumentError(
        "Constrained_Optimization.matroid_collinearity_system returned " *
        "unsupported value of type $(typeof(result))."
    ))
end


"""
    matroid_collinearity_system(n_points, collinear_sets; separating_bases=[])

Build Deepak's collinearity system with `hard_non_collinearity_relns=false`.
Optionally add one Rabinowitsch variable for each selected rank-three basis.
"""
function matroid_collinearity_system(
    n_points::Int,
    collinear_sets::Vector{NTuple{3,Int}};
    separating_bases::Vector{Vector{Int}}=Vector{Vector{Int}}(),
)
    raw_base = _matroid_collinearity_system(
        n_points,
        collinear_sets;
        hard_non_collinearity_relns=false,
    )

    base = _as_hc_system(raw_base)

    vars_xy = collect(variables(base))
    equations = Any[collect(expressions(base))...]

    if isempty(equations)
        push!(equations, 0.0 * vars_xy[1])
    end

    slack_variables = Variable[]
    for (index, basis) in enumerate(separating_bases)
        slack = Variable("u_sep_$index")
        push!(slack_variables, slack)
        determinant = _affine_minor_expression(vars_xy, basis)
        push!(equations, slack * determinant - 1.0)
    end

    return System(equations; variables=vcat(vars_xy, slack_variables))
end

function drawing_system(M; separating_bases=Vector{Vector{Int}}())
    Oscar.rank(M) == 3 || error("Polynomial drawing systems are used only for rank-three matroids.")
    return matroid_collinearity_system(
        length(M),
        collinearity_triples(M);
        separating_bases=separating_bases,
    )
end

function circle_guess(
    n_points::Int,
    bounds;
    phase::Real=0.0,
    jitter::Real=0.0,
    rng::Random.AbstractRNG=Random.default_rng(),
)
    (xmin, xmax), (ymin, ymax) = bounds
    radius = 0.35 * min(xmax - xmin, ymax - ymin)
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    guess = Float64[]
    for i in 1:n_points
        angle = phase + 2pi * (i - 1) / n_points
        x = center_x + radius * cos(angle)
        y = center_y + radius * sin(angle)

        if jitter > 0
            x += jitter * radius * randn(rng)
            y += jitter * radius * randn(rng)
        end

        push!(guess, x, y)
    end

    return guess
end

function append_rabinowitsch_coordinates(
    coordinates::AbstractVector,
    separating_bases::Vector{Vector{Int}})
    
    state = Float64.(collect(coordinates))

    for basis in separating_bases
        determinant = Float64(affine_minor_value(coordinates, basis))
        abs(determinant) > 1e-10 || throw(ArgumentError(
            "The initial affine guess has a nearly zero selected basis determinant for $basis."
        ))
        push!(state, inv(determinant))
    end

    return state
end

function rank_two_default_coordinates(n_points::Int, bounds)
    (xmin, xmax), (ymin, ymax) = bounds
    n_points >= 2 || error("A rank-two simple matroid needs at least two points.")

    xs = collect(range(xmin + 0.15 * (xmax - xmin), xmax - 0.15 * (xmax - xmin); length=n_points))
    y = (ymin + ymax) / 2

    coordinates = Float64[]
    for x in xs
        push!(coordinates, x, y)
    end
    return coordinates
end
