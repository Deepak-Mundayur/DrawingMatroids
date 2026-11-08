# Constrained_Optimization's repulsion between points and from the walls of `bounds`, plus a
# torque that opens up the angles between intersecting lines of the drawing.
function drawing_force_field(n_pts::Int, bounds, lines; k_p, k_wall, k_torque)
    repel = make_bounded_repelling_force(n_pts, bounds; k_p=k_p, k_wall=k_wall)
    return pos -> torque_force!(repel(pos), pos, lines, k_torque)
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