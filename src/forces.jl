# Constrained_Optimization's repulsion between points and from the walls of `bounds`, plus a
# torque that opens up the angles between intersecting lines of the drawing and a push of each
# point away from the lines it is not on. The push on each point is capped at `max_force`: the
# repulsion grows without bound as two points meet, and an uncapped push throws points out of
# the box, where the walls then hold them outside.
function drawing_force_field(n_pts::Int, bounds, lines; k_p, k_wall, k_torque, k_line, max_force)
    repel = make_bounded_repelling_force(n_pts, bounds; k_p=k_p, k_wall=k_wall)
    function field(pos)
        force = line_repulsion_force!(torque_force!(repel(pos), pos, lines, k_torque), pos, lines, n_pts, k_line)
        for p in 1:n_pts
            magnitude = hypot(force[2p - 1], force[2p])
            magnitude > max_force && (force[(2p - 1):(2p)] .*= max_force / magnitude)
        end
        return force
    end
    return field
end

# Push each point away from every drawn line (the segment between the line's outermost points)
# that it is not on, with strength k_line / distance, so no point seems to lie on such a line.
function line_repulsion_force!(force, pos, lines, n_pts, k_line)
    for flat in lines
        (ax, bx), (ay, by) = _line_segment(pos, flat)
        for p in 1:n_pts
            p in flat && continue
            px, py = pos[2p - 1], pos[2p]
            vx, vy = _segment_offset(px, py, ax, ay, bx, by)
            d2 = vx^2 + vy^2 + 1e-6
            force[2p - 1] += k_line * vx / d2
            force[2p] += k_line * vy / d2
        end
    end
    return force
end

# For every two lines meeting in exactly one point c, and every i on the first and j on the
# second line, push i and j perpendicular to their arms from c with strength k * cot(angle icj),
# the cotangent capped at max_cot, which turns the angle toward a right angle; c takes the
# reaction force.
function torque_force!(force, pos, lines, k_torque; max_cot=10.0)
    point(i) = (pos[2i - 1], pos[2i])
    add!(i, (fx, fy)) = (force[2i - 1] += fx; force[2i] += fy)
    for a in eachindex(lines), b in (a + 1):lastindex(lines)
        L1, L2 = lines[a], lines[b]
        shared = intersect(L1, L2)
        length(shared) == 1 || continue
        c = only(shared)
        for i in L1, j in L2
            (i == c || j == c) && continue
            (Ax, Ay), (Bx, By) = point(i) .- point(c), point(j) .- point(c)
            (Ax^2 + Ay^2 < 1e-6 || Bx^2 + By^2 < 1e-6) && continue
            cp = Ax * By - Ay * Bx
            cot = clamp((Ax * Bx + Ay * By) / (cp >= 0 ? max(cp, 1e-6) : min(cp, -1e-6)), -max_cot, max_cot)
            Fi, Fj = (-k_torque * cot) .* (-Ay, Ax), (k_torque * cot) .* (-By, Bx)
            add!(i, Fi)
            add!(j, Fj)
            add!(c, .-(Fi .+ Fj))
        end
    end
    return force
end
