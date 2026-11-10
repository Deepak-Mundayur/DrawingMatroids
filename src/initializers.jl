# n points evenly spaced on a circle in the box, rotated by `phase`; with `jitter`, each is
# moved by jitter × radius in a random direction.
function circle_guess(n_points::Int, bounds; phase::Real=0.0, jitter::Real=0.0)
    (xmin, xmax), (ymin, ymax) = bounds
    radius = 0.35 * min(xmax - xmin, ymax - ymin)
    angles = phase .+ 2pi .* (0:n_points-1) ./ n_points
    P = [(xmin + xmax) / 2 .+ radius .* cos.(angles)'; (ymin + ymax) / 2 .+ radius .* sin.(angles)']
    jitter > 0 && (P .+= jitter * radius .* [randn() for _ in 1:2, _ in 1:n_points])   # one draw per entry, x1, y1, x2, ...
    return vec(P)
end

# Points evenly spaced on a horizontal line (one point: the centre), for ranks 0 to 2.
function low_rank_default_coordinates(n_points::Int, bounds)
    (xmin, xmax), (ymin, ymax) = bounds
    n_points == 1 && return [(xmin + xmax) / 2, (ymin + ymax) / 2]
    xs = range(xmin + 0.15 * (xmax - xmin), xmax - 0.15 * (xmax - xmin); length=n_points)
    return [c for x in xs for c in (x, (ymin + ymax) / 2)]
end

# Affine drawing coordinates from a real sample of `C`.
# Component samples are validated only after this normalization: a raw sample is an
# arbitrary projective representative whose minors can be tiny purely from column scaling.
# A component of rank r inside a space of larger matrices (e.g. all points on a line in a
# rank-3 space) is sampled with extra rows; its column space gives an r×n matrix with the
# same column dependencies.
function start_coordinates(C::RealizationSpaces.RealizationComponent, reduction; bounds, real_tol)
    A = RealizationSpaces.sample(C; real_point=true)
    r = Oscar.rank(reduction.matroid)
    if size(A, 1) > r
        maximum(abs, imag.(A)) <= real_tol || error("The sample of the component is not numerically real.")
        A = svd(real.(A)).U[:, 1:r]' * real.(A)
    end
    return realization_matrix_to_coordinates(A, reduction; bounds=bounds, real_tol=real_tol)
end

# Cyclic flats of M grouped by rank 0..3: the input format of Maximat (Liwski, Mohammadi,
# Prébet, Def. 3.14), whose search runs in rank at most four.
function cyclic_flats_by_rank(M)
    Oscar.rank(M) <= 3 || error("The maximal-degeneration adapter is restricted to rank at most 3.")
    levels = [Set{Int}[] for _ in 1:4]
    for F in Oscar.cyclic_flats(M)
        isempty(F) || push!(levels[Oscar.rank(M, F) + 1], Set{Int}(F))
    end
    return levels
end

# Rank of X in the degeneration given by Maximat's reduced cyclic flats CF: drop loops,
# map elements to their parallel-class representatives, then take the valuation
# min(|Y|, 4, |Y \\ F| + rank(F)) over the type-2 and type-3 edges F (Def. 4.3).
function degeneration_rank(X, CF)
    loops = isempty(CF[1]) ? Set{Int}() : only(CF[1])
    representative = Dict(e => minimum(c) for c in CF[2] for e in c)
    Y = Set(get(representative, x, x) for x in X if !(x in loops))
    values = [length(Y), 4]
    for (r, level) in ((2, CF[3]), (3, CF[4])), F in level
        push!(values, r + length(setdiff(Y, F)))
    end
    return minimum(values)
end

"""
    maximal_degeneration_basis_transversal(M)

Maximal degenerations of `M` (rank 2 or 3) and a set of separating bases: bases of `M` such
that every maximal degeneration makes at least one of them dependent. Nonvanishing of their
determinants on the circuit variety therefore cuts out the realization space of `M`.
"""
function maximal_degeneration_basis_transversal(M)
    r = Oscar.rank(M)
    r in (2, 3) || error("Maximal-degeneration initialization supports ranks 2 and 3.")
    degenerations = Maximat.maximal_degenerations(cyclic_flats_by_rank(M), length(M); S=1:r, n_procs=Base.Threads.nthreads())

    bases = [Set(Int.(B)) for B in Oscar.bases(M)]
    killed = [Set(B for B in bases if degeneration_rank(B, N) < r) for N in degenerations]
    any(isempty, killed) && error("A maximal degeneration could not be separated from M by a basis of M.")
    separating_bases = [sort(collect(B)) for B in RealizationSpaces._greedy_transversal(killed)]
    return separating_bases, degenerations
end

# The reduced Rabinowitsch system of a simple rank-three matroid: its drawing system plus one
# Rabinowitsch variable per separating basis (instead of one per basis).
function rabinowitsch_system(simple)
    separating_bases, degenerations = maximal_degeneration_basis_transversal(simple)
    println("[draw_matroid]   $(length(degenerations)) maximal degenerations, $(length(separating_bases)) separating bases")
    return (system=drawing_system(simple; separating_bases), separating_bases)
end

# Project the circle guess onto the reduced Rabinowitsch system; return drawing coordinates
# that realize the matroid, or throw.
function project_rabinowitsch(rabinowitsch, coordinates, reduction; tol)
    minors = [det(affine_points_matrix(coordinates)[:, B]) for B in rabinowitsch.separating_bases]
    all(m -> abs(m) > 1e-10, minors) || error("The circle guess has a nearly zero separating basis determinant.")
    point, converged = project_onto_manifold(rabinowitsch.system, [coordinates; inv.(minors)]; tol=1e-10, max_iters=100)
    converged || error("The reduced Rabinowitsch projection did not converge.")
    coordinates = Float64.(point[1:length(coordinates)])
    validation = compute_drawing_validation(coordinates, reduction; tol)
    validation.valid || error("The reduced Rabinowitsch projection does not realize the matroid: $validation")
    return coordinates
end

# A real point of Oscar's exact realization space of M, as a Float64 matrix.
# - No variables: the realization is unique and rational.
# - No equations left: the variables are free, so random rationals at which no inequation
#   vanishes give an exact rational realization.
# - Dimension 0: the real solutions are isolated exactly; none means no real realization.
# Returns `nothing` when equations remain in positive dimension or no point was found, and a
# string saying why when M has no real realization. It never calls Oscar's `is_realizable`
# (a prime decomposition that can run for hours); it uses only what `realization_space`
# already found.
function exact_real_realization(M; tries::Int=100)
    X = Oscar.realization_space(M; char=0)
    Oscar.has_attribute(X, :is_realizable) && !Oscar.is_realizable(X) &&
        return "it is not realizable over any field of characteristic 0"
    A = Oscar.realization_matrix(X)
    at(p) = [Float64(Oscar.evaluate(A[i, j], p)) for i in 1:Oscar.nrows(A), j in 1:Oscar.ncols(A)]
    R = Oscar.ambient_ring(X)
    R isa Oscar.MPolyRing || return [Float64(Oscar.QQ(A[i, j])) for i in 1:Oscar.nrows(A), j in 1:Oscar.ncols(A)]
    I = Oscar.defining_ideal(X)
    if Oscar.is_zero(I)
        for _ in 1:tries
            p = [Oscar.QQ(rand(-30:30), rand(1:7)) for _ in 1:Oscar.ngens(R)]
            all(f -> !Oscar.is_zero(Oscar.evaluate(f, p)), Oscar.inequations(X)) || continue
            return at(p)
        end
        return nothing
    end
    Oscar.dim(I) == 0 || return nothing
    for s in first(Oscar.real_solutions(I))
        p = Float64.(s)
        all(f -> abs(Oscar.evaluate(f, p)) > 1e-8, Oscar.inequations(X)) || continue
        return at(p)
    end
    return "it has no real realization (Oscar's exact realization space is finite with no real point)"
end
