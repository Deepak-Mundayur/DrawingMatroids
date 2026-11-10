# Readability: a valid drawing is readable when its points are apart, off the lines they are not
# on, and its labels fit; the score, the projective view search, and the Ipopt polish.

# Readability thresholds, as fractions of the box width.
const _MIN_POINT_GAP = 0.10        # between two points
const _MIN_LINE_GAP = 0.06         # between a point and a drawn line it is not on

# Smallest distance from a point to a drawn line it is not on.
function _line_clearance(coords, reduction)
    clearance = Inf
    for flat in reduction.lines
        (ax, bx), (ay, by) = _line_segment(coords, flat)
        for p in eachindex(reduction.parallel_classes)
            p in flat || (clearance = min(clearance, hypot(_segment_offset(coords[2p - 1], coords[2p], ax, ay, bx, by)...)))
        end
    end
    return clearance
end

# Readability of a drawing inside `bounds`: the smallest of (closest two points) / _MIN_POINT_GAP
# and (closest point to a drawn line it is not on) / _MIN_LINE_GAP, in box widths; 0 if a label
# leaves the box or overlaps a point or another label. A drawing is readable when this is >= 1.
function _readability(coords, reduction, bounds)
    width = bounds[1][2] - bounds[1][1]
    k = length(reduction.parallel_classes)
    point(i) = (coords[2i - 1], coords[2i])
    score = Inf
    for i in 1:k, j in (i + 1):k
        score = min(score, hypot((point(i) .- point(j))...) / (_MIN_POINT_GAP * width))
    end
    score = min(score, _line_clearance(coords, reduction) / (_MIN_LINE_GAP * width))
    boxes = [spot[5] for spot in values(_label_layout(coords, reduction, bounds))]
    marker = 0.015 * width
    for (a, b) in enumerate(boxes)
        _box_outside(b, bounds) > 0 && return 0.0
        any(_box_overlap(b, q) > 0 for q in boxes[(a + 1):end]) && return 0.0
    end
    # A label's own point sits 0.03 widths from its box; no point may come within a marker radius.
    any(_box_distance(b, point(i)) < marker for b in boxes for i in 1:k) && return 0.0
    return score
end

# A projective transformation keeps every line and every non-collinearity, so it turns a
# drawing into another drawing of the same matroid; it acts by an invertible matrix on the left,
# under which every component of a realization space is invariant, so it also stays on the
# same component. A readable drawing is kept; otherwise the frame with the largest `_readability` is
# chosen (random search, then local refinement) and used if it still validates.
function _readable_frame(coords, reduction, bounds; tol, candidates=200, refinements=400)
    k = length(reduction.parallel_classes)
    k >= 2 || return coords
    H = [ones(1, k); reshape(coords[1:2k], 2, k)]
    function framed(G)
        B = G * H
        all(abs.(B[1, :]) .> 1e-6) || return nothing
        return _fit_coordinates_to_bounds(vec(B[2:3, :] ./ B[1:1, :]), bounds)
    end
    score(c) = isnothing(c) ? -Inf : _readability(c, reduction, bounds)
    original = score(framed(Matrix(1.0I, 3, 3)))
    original >= 1 && return coords
    best_G, best = Matrix(1.0I, 3, 3), original
    for _ in 1:candidates
        G = randn(3, 3)
        s = score(framed(G))
        s > best && ((best_G, best) = (G, s))
    end
    step = 0.3
    for _ in 1:refinements
        best >= 1.5 && break
        G = best_G .+ step * sqrt(sum(abs2, best_G) / 9) .* randn(3, 3)
        s = score(framed(G))
        s > best ? ((best_G, best) = (G, s)) : (step *= 0.995)
    end
    best > original || return coords
    c = framed(best_G)
    return compute_drawing_validation(c, reduction; tol).valid ? c : coords
end

# Readability as a nonlinear program, solved by Ipopt from a valid drawing. Unknowns: the
# points in [0,1]^2 and s. Maximize s subject to
#   every nonbasis determinant = 0                       (lines stay lines)
#   |p_i - p_j|^2   >= s * _MIN_POINT_GAP^2              (points apart)
#   dist(p, line)^2 >= s * _MIN_LINE_GAP^2, p not on it  (points off lines they are not on)
# so sqrt(s) is the readability without the label rule. The distance is to the whole line,
# which is stricter than to the drawn segment.
function _ipopt_layout(coords, reduction; max_iter=3000, time_limit=60.0)
    k = length(reduction.parallel_classes)
    start = _fit_coordinates_to_bounds(Float64.(coords[1:2k]), ((0.0, 1.0), (0.0, 1.0)))
    model = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_attribute(model, "sb", "yes")
    JuMP.set_attribute(model, "max_iter", max_iter)
    JuMP.set_attribute(model, "max_cpu_time", time_limit)
    JuMP.set_attribute(model, "tol", 1e-10)
    JuMP.set_attribute(model, "constr_viol_tol", 1e-11)
    JuMP.@variable(model, 0 <= x[i=1:k] <= 1, start = start[2i - 1])
    JuMP.@variable(model, 0 <= y[i=1:k] <= 1, start = start[2i])
    det3(a, b, c) = (x[b] - x[a]) * (y[c] - y[a]) - (x[c] - x[a]) * (y[b] - y[a])
    for NB in Oscar.nonbases(reduction.simple_matroid)
        a, b, c = sort(Int.(collect(NB)))
        JuMP.@constraint(model, det3(a, b, c) == 0)
    end
    gap = minimum(((start[2i - 1] - start[2j - 1])^2 + (start[2i] - start[2j])^2) / _MIN_POINT_GAP^2 for i in 1:k for j in (i + 1):k)
    JuMP.@variable(model, 0 <= s <= 25, start = min(gap, 1.0))
    for i in 1:k, j in (i + 1):k
        JuMP.@constraint(model, (x[i] - x[j])^2 + (y[i] - y[j])^2 >= s * _MIN_POINT_GAP^2)
    end
    for L in reduction.lines
        a, b = L[1], L[2]
        for p in 1:k
            p in L || JuMP.@constraint(model, det3(a, b, p)^2 >= s * _MIN_LINE_GAP^2 * ((x[a] - x[b])^2 + (y[a] - y[b])^2))
        end
    end
    JuMP.@objective(model, Max, s)
    JuMP.optimize!(model)
    return [c for i in 1:k for c in (JuMP.value(x[i]), JuMP.value(y[i]))]
end

# The polished drawing, kept only if it still realizes the matroid and is more readable.
function _ipopt_polish(coords, reduction, bounds; tol)
    length(reduction.parallel_classes) >= 3 || return coords
    k = length(reduction.parallel_classes)
    before = _readability(_fit_coordinates_to_bounds(Float64.(coords[1:2k]), bounds), reduction, bounds)
    polished = try
        _fit_coordinates_to_bounds(_ipopt_layout(coords, reduction), bounds)
    catch err
        err isa InterruptException && rethrow()
        nothing
    end
    after = isnothing(polished) || !compute_drawing_validation(polished, reduction; tol).valid ? -Inf :
            _readability(polished, reduction, bounds)
    return after > before ? polished : coords
end
