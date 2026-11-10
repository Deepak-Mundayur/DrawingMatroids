# Remove loops and merge parallel classes. The simple matroid has one element per class,
# labelled 1:k in the order of the classes' smallest elements; `lines` are its drawn lines.
function drawing_reduction(M)
    n = length(M)
    sort(Int.(collect(Oscar.matroid_groundset(M)))) == 1:n ||
        throw(ArgumentError("draw_matroid requires the matroid's ground set to be 1:$n."))
    loops = Set(Int.(Oscar.loops(M)))
    classes = unique([sort(setdiff(Int.(Oscar.closure(M, [i])), loops)) for i in 1:length(M) if !(i in loops)])
    representatives = first.(classes)
    index = Dict(e => j for (j, e) in enumerate(representatives))
    restricted = Oscar.restriction(M, representatives)
    simple = Oscar.matroid_from_bases([[index[e] for e in B] for B in Oscar.bases(restricted)], length(classes))
    return (matroid=M, simple_matroid=simple, parallel_classes=classes, lines=rank_two_flats(simple))
end

# The r×n realization matrix of a drawing: column [1, x, y] (its first r entries) for every
# element of a parallel class, zero for loops.
function coordinates_to_matrix(coords::AbstractVector, reduction)
    M = reduction.matroid
    k = length(reduction.parallel_classes)
    length(coords) >= 2k || throw(ArgumentError("Expected at least $(2k) affine coordinates, got $(length(coords))."))
    A = zeros(Oscar.rank(M), length(M))
    for (c, class) in enumerate(reduction.parallel_classes)
        A[:, class] .= [1.0, real(coords[2c - 1]), real(coords[2c])][1:size(A, 1)]
    end
    return A
end

# Centre the points (x1, y1, x2, ...) in the box `bounds` and scale them, keeping the aspect
# ratio, so they span 70% of it.
function _fit_coordinates_to_bounds(coords::Vector{Float64}, bounds)
    isempty(coords) && return coords
    P = reshape(coords, 2, :)
    lo, hi = minimum(P; dims=2), maximum(P; dims=2)
    width = [b[2] - b[1] for b in bounds]
    scale = minimum(s > 1e-12 ? 0.7w / s : 1.0 for (s, w) in zip(hi .- lo, width))
    center = [(b[1] + b[2]) / 2 for b in bounds]
    return vec(center .+ scale .* (P .- (lo .+ hi) ./ 2))
end

# Affine drawing coordinates of a real r×n realization matrix: one point per parallel class,
# in an affine chart (the identity, else a random invertible G with every point off the line
# at infinity), fitted to `bounds`.
function realization_matrix_to_coordinates(A::AbstractMatrix, reduction; bounds, real_tol::Real=1e-7)
    r, n = Oscar.rank(reduction.matroid), length(reduction.matroid)
    size(A) == (r, n) || throw(ArgumentError("Realization matrix has size $(size(A)); expected ($r, $n)."))
    maximum(abs, imag.(A); init=0.0) <= real_tol || error("The supplied realization matrix is not numerically real.")
    B = real.(A)[:, first.(reduction.parallel_classes)]
    isempty(B) && return Float64[]
    for attempt in 0:100
        G = attempt == 0 ? Matrix(1.0I, r, r) : randn(r, r)
        C = G * B
        abs(det(G)) > 1e-9 && all(abs.(C[1, :]) .> 1e-9) || continue
        return _fit_coordinates_to_bounds(vec([C[2:r, :] ./ C[1:1, :]; zeros(3 - r, size(C, 2))]), bounds)
    end
    error("Could not find a real affine chart containing every nonloop point.")
end

# Lines of the drawing: the whole ground set in rank 2, the rank-2 cyclic flats in rank 3.
function rank_two_flats(M)
    Oscar.rank(M) == 2 && return [collect(1:length(M))]
    Oscar.rank(M) == 3 || return Vector{Vector{Int}}()
    return [sort(Int.(F)) for F in Oscar.cyclic_flats(M, 2)]
end
