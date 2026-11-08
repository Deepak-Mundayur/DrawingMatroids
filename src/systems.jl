# Symbolic or numeric 3×n matrix [1; x; y] of affine points, coords = (x1, y1, x2, y2, ...).
affine_points_matrix(coords::AbstractVector) = [ones(Int, 1, length(coords) ÷ 2); reshape(coords, 2, :)]

"""
    drawing_system(M; separating_bases=[])

Affine collinearity system of a simple rank-three matroid: one minor per nonbasis,
plus one Rabinowitsch equation `u * det(B) - 1` per separating basis `B`.
"""
function drawing_system(M; separating_bases=Vector{Vector{Int}}())
    Oscar.rank(M) == 3 || error("Polynomial drawing systems are used only for rank-three matroids.")
    @var p[1:2, 1:length(M)] u[1:length(separating_bases)]
    A = affine_points_matrix(vec(p))
    equations = Expression[RealizationSpaces._nonbasis_equations(A, M);
                           [u[i] * det(A[:, B]) - 1 for (i, B) in enumerate(separating_bases)]]
    isempty(equations) && push!(equations, 0 * p[1])
    return System(equations; variables=[vec(p); u])
end

function circle_guess(
    n_points::Int,
    bounds;
    phase::Real=0.0,
    jitter::Real=0.0,
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
            x += jitter * radius * randn()
            y += jitter * radius * randn()
        end

        push!(guess, x, y)
    end

    return guess
end

function rank_two_default_coordinates(n_points::Int, bounds)
    (xmin, xmax), (ymin, ymax) = bounds
    n_points >= 2 || error("A rank-two simple matroid needs at least two points.")
    xs = range(xmin + 0.15 * (xmax - xmin), xmax - 0.15 * (xmax - xmin); length=n_points)
    return [c for x in xs for c in (x, (ymin + ymax) / 2)]
end
