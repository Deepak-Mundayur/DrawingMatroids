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
