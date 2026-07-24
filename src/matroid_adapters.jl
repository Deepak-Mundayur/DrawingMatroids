const _DEFAULT_CHART_TOL = 1e-9
const _DEFAULT_CHART_ATTEMPTS = 100

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

function _parallel_classes(M)
    n = length(M)
    loops = Set(Int.(collect(Oscar.loops(M))))
    seen = Set{Int}()
    classes = Vector{Vector{Int}}()

    for i in 1:n
        (i in loops || i in seen) && continue

        cls = Int[]
        for j in 1:n
            j in loops && continue
            if Oscar.rank(M, [i, j]) == 1
                push!(cls, j)
            end
        end

        sort!(cls)
        push!(classes, cls)
        union!(seen, cls)
    end

    sort!(classes; by=first)
    return classes, sort!(collect(loops))
end

function _simple_matroid_from_classes(M, classes::Vector{Vector{Int}})
    r = Oscar.rank(M)
    r == 0 && return Oscar.uniform_matroid(0, length(classes))

    class_of = Dict{Int,Int}()

    for (i, cls) in enumerate(classes), element in cls
        class_of[element] = i
    end

    basis_keys = Set{Tuple}()
    for B in Oscar.bases(M)
        mapped = sort(unique([class_of[Int(e)] for e in collect(B)]))
        length(mapped) == r && push!(basis_keys, Tuple(mapped))
    end

    simple_bases = [collect(key) for key in basis_keys]
    isempty(simple_bases) && error("Could not construct the simplified matroid from its bases.")

    return Oscar.matroid_from_bases(simple_bases, length(classes))
end

function drawing_reduction(M)
    _require_groundset_1_to_n(M)

    classes, loops = _parallel_classes(M)
    simple = _simple_matroid_from_classes(M, classes)
    representatives = first.(classes)
    grouped_labels = [Any[e for e in cls] for cls in classes]
    loop_labels = Any[e for e in loops]

    return DrawingReduction(M, simple, classes, representatives, grouped_labels, loops, loop_labels)
end

function coordinate_count(reduction::DrawingReduction)
    return 2 * length(reduction.parallel_classes)
end

function coordinates_to_matrix(coords::AbstractVector, reduction::DrawingReduction)
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

        column = if r == 3
            [1.0, x, y]
        elseif r == 2
            [1.0, x]
        else
            error("Drawing coordinates are implemented only for ranks 2 and 3.")
        end

        for element in cls
            A[:, element] .= column
        end
    end

    # Loop columns remain zero.
    return A
end

function _random_invertible_matrix(r::Int, rng::Random.AbstractRNG)
    for _ in 1:100
        G = randn(rng, r, r)
        abs(det(G)) > _DEFAULT_CHART_TOL &&
            return G
    end
    error("Could not generate an invertible chart matrix.")
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
    A::AbstractMatrix, reduction::DrawingReduction;
    bounds=((-5.0, 5.0), (-5.0, 5.0)), real_tol::Real=1e-7,
    rng::Random.AbstractRNG=Random.default_rng(),
)
    M = reduction.matroid
    r = Oscar.rank(M)
    n = length(M)

    size(A) == (r, n) ||
        throw(ArgumentError("Realization matrix has size $(size(A)); expected ($r, $n)."))

    maximum(abs.(imag.(ComplexF64.(A)))) <= real_tol ||
        error("The supplied realization matrix is not numerically real.")

    Ar = Matrix{Float64}(real.(A))
    reps = reduction.representatives

    candidates = Matrix{Float64}[Matrix{Float64}(I, r, r)]
    for _ in 1:_DEFAULT_CHART_ATTEMPTS
        push!(candidates, _random_invertible_matrix(r, rng))
    end

    for G in candidates
        B = G * Ar
        all(abs(B[1, j]) > _DEFAULT_CHART_TOL for j in reps) || continue

        coords = Float64[]
        for j in reps
            normalized = B[:, j] ./ B[1, j]
            if r == 3
                push!(coords, normalized[2], normalized[3])
            elseif r == 2
                push!(coords, normalized[2], 0.0)
            else
                error("Drawing is implemented only for ranks 2 and 3.")
            end
        end

        return _fit_coordinates_to_bounds(coords, bounds)
    end

    error(
        "Could not find a real affine chart containing every nonloop point."
    )
end

function _matrix_from_user_input(data, reduction::DrawingReduction)
    M = reduction.matroid
    r = Oscar.rank(M)
    n = length(M)
    k = length(reduction.parallel_classes)

    if data isa AbstractMatrix && size(data) == (r, n)
        return Matrix(data)
    elseif data isa AbstractMatrix && size(data) == (r, k)
        simple_matrix = Matrix(data)
        expanded = zeros(eltype(simple_matrix), r, n)
        for (class_index, cls) in enumerate(reduction.parallel_classes)
            for element in cls
                expanded[:, element] .= simple_matrix[:, class_index]
            end
        end
        return expanded
    end

    return nothing
end

function _require_numerically_real(data; tol::Real=1e-7)
    values = ComplexF64.(collect(data))
    isempty(values) && return nothing
    maximum(abs.(imag.(values))) <= tol || throw(ArgumentError(
        "Affine drawing coordinates must be numerically real."
    ))
    return nothing
end

function user_input_to_coordinates(
    data,
    reduction::DrawingReduction;
    bounds=((-5.0, 5.0), (-5.0, 5.0)),
    real_tol::Real=1e-7,
    rng::Random.AbstractRNG=Random.default_rng(),
)
    data === nothing && throw(ArgumentError(
        "start_point_strategy=:user requires `start_point` to be a realization matrix or affine coordinates."
    ))

    matrix_input = _matrix_from_user_input(data, reduction)
    if matrix_input !== nothing
        return realization_matrix_to_coordinates(
            matrix_input,
            reduction;
            bounds=bounds,
            real_tol=real_tol,
            rng=rng,
        )
    end

    k = length(reduction.parallel_classes)
    n = length(reduction.matroid)

    if data isa AbstractMatrix
        _require_numerically_real(data; tol=real_tol)
        B = Matrix(data)

        if size(B) == (2, k)
            return _fit_coordinates_to_bounds(vec(B), bounds)
        elseif size(B) == (k, 2)
            return _fit_coordinates_to_bounds(vec(permutedims(B)), bounds)
        elseif size(B) == (2, n)
            coords = Float64[]
            for rep in reduction.representatives
                push!(coords, Float64(real(B[1, rep])), Float64(real(B[2, rep])))
            end
            return _fit_coordinates_to_bounds(coords, bounds)
        elseif size(B) == (n, 2)
            coords = Float64[]
            for rep in reduction.representatives
                push!(coords, Float64(real(B[rep, 1])), Float64(real(B[rep, 2])))
            end
            return _fit_coordinates_to_bounds(coords, bounds)
        end
    elseif data isa AbstractVector
        _require_numerically_real(data; tol=real_tol)
        values = Float64.(real.(collect(data)))
        if length(values) == 2 * k
            return _fit_coordinates_to_bounds(values, bounds)
        elseif length(values) == 2 * n
            coords = Float64[]
            for rep in reduction.representatives
                push!(coords, values[2 * rep - 1], values[2 * rep])
            end
            return _fit_coordinates_to_bounds(coords, bounds)
        end
    end

    throw(ArgumentError(
        "Unsupported user start point. Supply an r×n realization matrix, a 2×k/k×2 affine matrix, or a length-2k affine vector."
    ))
end

function _push_unique_set!(sets::Vector{Set{Int}}, candidate::Set{Int})
    any(existing == candidate for existing in sets) || push!(sets, candidate)
    return sets
end

"""
    cyclic_flats_by_rank(M)

Return the four rank levels consumed by Rémi's translated hypergraph code. The
actual cyclic flats are included, and the full ground set is also inserted at
`rank(M)` as a redundant total-rank constraint. That extra edge is necessary
for rank-two/rank-three matroids with coloops because the original program was
written with the ambient rank fixed externally.
"""
function cyclic_flats_by_rank(M)
    Oscar.rank(M) <= 3 || error("Rémi's drawing adapter is restricted to rank at most 3.")

    levels = [Set{Int}[] for _ in 1:4]
    loops = Set(Int.(collect(Oscar.loops(M))))
    !isempty(loops) && push!(levels[1], loops)

    for C in Oscar.circuits(M)
        circuit = Int.(collect(C))
        flat = Set(Int.(collect(Oscar.closure(M, circuit))))
        flat_rank = Oscar.rank(M, collect(flat))
        0 <= flat_rank <= 3 && _push_unique_set!(levels[flat_rank + 1], flat)
    end

    # Rémi's rank-four implementation treats the ambient rank as fixed outside
    # the cyclic-flat data. For ranks two and three we insert the redundant
    # total-rank constraint E at its actual rank, even when E is not cyclic
    # because the matroid has coloops. This prevents the translated hypergraph
    # search from admitting higher-rank candidates.
    E = Set(1:length(M))
    _push_unique_set!(levels[Oscar.rank(M) + 1], E)

    for level in levels
        sort!(level; by=s -> (length(s), Tuple(sort(collect(s)))))
    end

    return levels
end

function cyclic_total_rank(CF, original_rank::Int, d::Int)
    E = Set(1:d)
    total_rank = original_rank

    for rank_value in 0:3
        if any(F -> F == E, CF[rank_value + 1])
            total_rank = min(total_rank, rank_value)
        end
    end

    return total_rank
end

function rank_from_cyclic_flats(X, CF, original_rank::Int, d::Int)
    subset = Set(Int.(collect(X)))
    values = Int[length(subset)] # the empty cyclic flat

    for rank_value in 0:3, F in CF[rank_value + 1]
        push!(values, rank_value + length(setdiff(subset, F)))
    end

    total_rank = cyclic_total_rank(CF, original_rank, d)
    push!(values, total_rank)
    return minimum(values)
end

function maximal_degeneration_basis_transversal(
    M;
    preprocess::Bool=true,
    verbosity::Int=0,
    n_procs::Int=Base.Threads.nthreads(),
)
    r = Oscar.rank(M)
    d = length(M)
    r in (2, 3) || error("Maximal-degeneration initialization supports ranks 2 and 3.")

    CF = cyclic_flats_by_rank(M)
    degenerations = Maximat.maximal_degenerations(
        CF,
        d;
        S=collect(1:r),
        v=verbosity,
        preprocess=preprocess,
        n_procs=n_procs,
    )

    target_bases = sort(
        [sort(Int.(collect(B))) for B in Oscar.bases(M)];
        by=B -> Tuple(B),
    )

    degeneration_basis_sets = Vector{Set{Tuple}}()
    for degeneration in degenerations
        killed = Set{Tuple}()
        for B in target_bases
            if rank_from_cyclic_flats(B, degeneration, r, d) < r
                push!(killed, Tuple(B))
            end
        end
        isempty(killed) && error(
            "A maximal degeneration could not be separated from the target by a target basis."
        )
        push!(degeneration_basis_sets, killed)
    end

    remaining = copy(degeneration_basis_sets)
    selected = Vector{Vector{Int}}()

    while !isempty(remaining)
        counts = Dict{Tuple,Int}()
        for killed in remaining, basis in killed
            counts[basis] = get(counts, basis, 0) + 1
        end

        keys_sorted = sort(collect(keys(counts)); by=identity)
        best = first(keys_sorted)
        best_count = counts[best]

        for candidate in keys_sorted
            candidate_count = counts[candidate]
            if candidate_count > best_count
                best = candidate
                best_count = candidate_count
            end
        end

        push!(selected, collect(best))
        filter!(killed -> !(best in killed), remaining)
    end

    return selected, degenerations
end

function rank_two_flats(M)
    r = Oscar.rank(M)
    n = length(M)

    if r == 2
        return n >= 2 ? [collect(1:n)] : Vector{Vector{Int}}()
    elseif r != 3
        return Vector{Vector{Int}}()
    end

    flats = Set{Int}[]
    for NB in Oscar.nonbases(M)
        triple = sort(Int.(collect(NB)))
        length(triple) == 3 || continue
        flat = Set(Int.(collect(Oscar.closure(M, triple))))
        length(flat) >= 3 && _push_unique_set!(flats, flat)
    end

    return [sort(collect(F)) for F in flats]
end
