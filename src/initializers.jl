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
function maximal_degeneration_basis_transversal(M; n_procs::Int=Base.Threads.nthreads())
    r = Oscar.rank(M)
    r in (2, 3) || error("Maximal-degeneration initialization supports ranks 2 and 3.")
    degenerations = Maximat.maximal_degenerations(cyclic_flats_by_rank(M), length(M); S=1:r, n_procs=n_procs)

    bases = [Set(Int.(B)) for B in Oscar.bases(M)]
    killed = [Set(B for B in bases if degeneration_rank(B, N) < r) for N in degenerations]
    any(isempty, killed) && error("A maximal degeneration could not be separated from M by a basis of M.")
    separating_bases = [sort(collect(B)) for B in RealizationSpaces._greedy_transversal(killed)]
    return separating_bases, degenerations
end

"""
    maximal_degeneration_start(reduction; bounds, attempts, tol)

Exact real start point for a simple rank-three matroid: compute its maximal
degenerations, choose a greedy transversal of separating bases, and project a circle
guess onto the reduced Rabinowitsch system (one Rabinowitsch variable per separating
basis instead of one per basis).
"""
function maximal_degeneration_start(reduction; bounds, attempts::Int, tol::Real)
    simple = reduction.simple_matroid
    point_count = length(simple)

    println("[DrawingMatroids] Computing maximal degenerations with $(Base.Threads.nthreads()) Julia thread(s)...")
    separating_bases, degenerations = maximal_degeneration_basis_transversal(simple)
    println("[DrawingMatroids] Found $(length(degenerations)) maximal degenerations and $(length(separating_bases)) separating bases.")

    system = drawing_system(simple; separating_bases=separating_bases)
    last_error = nothing

    for attempt in 1:attempts
        phase = attempt == 1 ? 0.0 : 2pi * rand()
        jitter = 0.015 * (attempt - 1)

        try
            coordinates = circle_guess(point_count, bounds; phase=phase, jitter=jitter)
            A = affine_points_matrix(coordinates)
            minors = [det(A[:, B]) for B in separating_bases]
            all(m -> abs(m) > 1e-10, minors) || error("The circle guess has a nearly zero separating basis determinant.")
            guess = [coordinates; inv.(minors)]

            println("[DrawingMatroids] Projecting onto the reduced Rabinowitsch system " *
                    "(attempt $attempt/$attempts, $(length(separating_bases)) separating bases)...")

            point, converged = project_onto_manifold(system, guess; tol=1e-10, max_iters=100)
            if !converged
                last_error = ErrorException("Reduced Rabinowitsch projection did not converge on attempt $attempt.")
                continue
            end

            coordinates = Float64.(point[1:(2 * point_count)])
            validation = compute_drawing_validation(coordinates, reduction; tol=tol)
            validation.valid && return coordinates
            last_error = ErrorException("Projected start point does not realize the matroid: $validation")
        catch err
            err isa InterruptException && rethrow()
            last_error = err
        end
    end

    error("Could not obtain an exact real start point from the reduced Rabinowitsch system " *
          "after $attempts attempts. Last error: $(sprint(showerror, last_error))")
end
