"""
    _require_groundset_1_to_n(M)

Throw an `ArgumentError` unless `M`'s ground set is exactly `1:length(M)`.
"""
function _require_groundset_1_to_n(M)
    n = length(M)
    groundset = sort(Int.(collect(Oscar.matroid_groundset(M))))
    groundset == collect(1:n) || throw(ArgumentError(
        "DrawingMatroids requires the matroid's ground set to be 1:$n; got $groundset."
    ))
    return nothing
end

# Remove loops and merge parallel classes. The simple matroid has one element per class,
# labelled 1:k in the order of the classes' smallest elements.
function drawing_reduction(M)
    _require_groundset_1_to_n(M)
    loops = Set(Int.(Oscar.loops(M)))
    classes = unique([sort(setdiff(Int.(Oscar.closure(M, [i])), loops)) for i in 1:length(M) if !(i in loops)])
    representatives = first.(classes)
    index = Dict(e => j for (j, e) in enumerate(representatives))
    restricted = Oscar.restriction(M, representatives)
    simple = Oscar.matroid_from_bases([[index[e] for e in B] for B in Oscar.bases(restricted)], length(classes))
    return (matroid=M, simple_matroid=simple, parallel_classes=classes)
end

function coordinates_to_matrix(coords::AbstractVector, reduction)
    M = reduction.matroid
    r = Oscar.rank(M)
    n = length(M)
    k = length(reduction.parallel_classes)

    length(coords) >= 2 * k ||
        throw(ArgumentError("Expected at least $(2 * k) affine coordinates, got $(length(coords))."))

    A = zeros(Float64, r, n)

    for (class_index, cls) in enumerate(reduction.parallel_classes)
        x = Float64(real(coords[2 * class_index - 1]))
        y = Float64(real(coords[2 * class_index]))

        column = [1.0, x, y][1:r]

        for element in cls
            A[:, element] .= column
        end
    end

    # Loop columns remain zero.
    return A
end

function _fit_coordinates_to_bounds(coords::Vector{Float64}, bounds)
    (xmin, xmax), (ymin, ymax) = bounds
    k = length(coords) ÷ 2
    k == 0 && return coords

    xs = [coords[2 * i - 1] for i in 1:k]
    ys = [coords[2 * i] for i in 1:k]

    cx = (minimum(xs) + maximum(xs)) / 2
    cy = (minimum(ys) + maximum(ys)) / 2
    span_x = maximum(xs) - minimum(xs)
    span_y = maximum(ys) - minimum(ys)

    available_x = 0.70 * (xmax - xmin)
    available_y = 0.70 * (ymax - ymin)
    scale_x = span_x > 1e-12 ? available_x / span_x : 1.0
    scale_y = span_y > 1e-12 ? available_y / span_y : 1.0
    scale = min(scale_x, scale_y)

    target_cx = (xmin + xmax) / 2
    target_cy = (ymin + ymax) / 2

    fitted = similar(coords)
    for i in 1:k
        fitted[2 * i - 1] = target_cx + scale * (xs[i] - cx)
        fitted[2 * i] = target_cy + scale * (ys[i] - cy)
    end
    return fitted
end

function realization_matrix_to_coordinates(
    A::AbstractMatrix, reduction;
    bounds, real_tol::Real=1e-7
)
    M = reduction.matroid
    r = Oscar.rank(M)
    n = length(M)

    size(A) == (r, n) ||
        throw(ArgumentError("Realization matrix has size $(size(A)); expected ($r, $n)."))

    maximum(abs.(imag.(ComplexF64.(A)))) <= real_tol ||
        error("The supplied realization matrix is not numerically real.")

    Ar = Matrix{Float64}(real.(A))
    reps = first.(reduction.parallel_classes)

    # An affine chart: the identity, else a random invertible G with every point off the line at infinity.
    for attempt in 0:100
        G = attempt == 0 ? Matrix(1.0I, r, r) : randn(r, r)
        abs(det(G)) > 1e-9 || continue
        B = G * Ar
        all(abs(B[1, j]) > 1e-9 for j in reps) || continue

        coords = Float64[]
        for j in reps
            normalized = B[:, j] ./ B[1, j]
            push!(coords, r >= 2 ? normalized[2] : 0.0, r == 3 ? normalized[3] : 0.0)
        end

        return _fit_coordinates_to_bounds(coords, bounds)
    end

    error(
        "Could not find a real affine chart containing every nonloop point."
    )
end

# Lines of the drawing: the whole ground set in rank 2, the rank-2 cyclic flats in rank 3.
function rank_two_flats(M)
    Oscar.rank(M) == 2 && return [collect(1:length(M))]
    Oscar.rank(M) == 3 || return Vector{Vector{Int}}()
    return [sort(Int.(F)) for F in Oscar.cyclic_flats(M, 2)]
end
