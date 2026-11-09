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
# which turns the angle icj toward a right angle; c takes the reaction force.
function torque_force!(force, pos, lines, k_torque; max_cot=10.0)
    n_lines = length(lines)
    
    for L1_idx in 1:n_lines
        for L2_idx in (L1_idx+1):n_lines
            L1 = lines[L1_idx]
            L2 = lines[L2_idx]
            
            shared_pts = intersect(L1, L2)
            if length(shared_pts) == 1
                c = shared_pts[1]
                idx_cx, idx_cy = 2*c - 1, 2*c
                cx, cy = pos[idx_cx], pos[idx_cy]
                
                for i in L1
                    if i == c continue end
                    idx_ix, idx_iy = 2*i - 1, 2*i
                    Ax, Ay = pos[idx_ix] - cx, pos[idx_iy] - cy
                    norm2_A = Ax^2 + Ay^2
                    if norm2_A < 1e-6 continue end 
                    
                    for j in L2
                        if j == c continue end
                        idx_jx, idx_jy = 2*j - 1, 2*j
                        Bx, By = pos[idx_jx] - cx, pos[idx_jy] - cy
                        norm2_B = Bx^2 + By^2
                        if norm2_B < 1e-6 continue end 
                        
                        dp = Ax * Bx + Ay * By
                        cp = Ax * By - Ay * Bx
                        
                        cp_safe = cp >= 0.0 ? max(cp, 1e-6) : min(cp, -1e-6)
                        
                        # Clamp the cotangent to prevent gradient explosion
                        cot_theta = clamp(dp / cp_safe, -max_cot, max_cot)
                        
                        # Tangential force on point i
                        fti = -k_torque * cot_theta 
                        Fix = fti * (-Ay) 
                        Fiy = fti * (Ax)
                        force[idx_ix] += Fix
                        force[idx_iy] += Fiy
                        
                        # Equal and opposite tangential force on point j
                        ftj = k_torque * cot_theta
                        Fjx = ftj * (-By)
                        Fjy = ftj * (Bx)
                        force[idx_jx] += Fjx
                        force[idx_jy] += Fjy
                        
                        # Apply reaction force to pivot 'c' to conserve momentum
                        force[idx_cx] -= (Fix + Fjx)
                        force[idx_cy] -= (Fiy + Fjy)
                    end
                end
            end
        end
    end
    return force
end