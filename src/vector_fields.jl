function make_bounded_repelling_force(n_pts::Int, bounds; k_p=10.0, k_wall=1000000000000000.0)
    (xmin, xmax), (ymin, ymax) = bounds

    function V(pos)
        force = zeros(length(pos))
        
        for i in 1:n_pts
            idx_x = 2*i - 1
            idx_y = 2*i
            x, y = pos[idx_x], pos[idx_y]
            
            #  Point-to-Point Repulsion 
            for j in 1:n_pts
                if i != j
                    jx, jy = pos[2*j - 1], pos[2*j]
                    dx = x - jx
                    dy = y - jy
                    dist_sq = dx^2 + dy^2 + 1e-6
                    
                    force[idx_x] += k_p * dx / dist_sq
                    force[idx_y] += k_p * dy / dist_sq
                end
            end
            
            # The bounding box forces ( Inverse square barrier)
            force[idx_x] += k_wall / (x - xmin)^2
            force[idx_x] -= k_wall / (xmax - x)^2
            
            force[idx_y] += k_wall / (y - ymin)^2
            force[idx_y] -= k_wall / (ymax - y)^2
        end
        
        # If the system has a saturation slack variable at the end, we assume that the its component of the vector field in that direction is 0 
        return force
    end
    
    return V
end



function make_bounded_repelling_force(n_pts::Int, bounds, collinearity_tuples::Vector{NTuple{3, Int}}; k_p=10.0, k_wall=100.0, k_torque=0.5)
    
    collinearity_relns = collect.(collinearity_tuples)
    
    function V(pos)
        force = zeros(length(pos))
        
        point_repulsion_force!(force, pos, n_pts, k_p)
        bounding_box_force!(force, pos, n_pts, bounds, k_wall)
        torque_force!(force, pos, collinearity_relns, k_torque)
        
        # If the system has a saturation slack variable at the end, 
        # we assume that its component of the vector field in that direction is 0 
        return force
    end
    
    return V
end

function point_repulsion_force!(force, pos, n_pts, k_p)
    for i in 1:n_pts
        idx_x, idx_y = 2*i - 1, 2*i
        x, y = pos[idx_x], pos[idx_y]
        
        for j in 1:n_pts
            if i != j
                jx, jy = pos[2*j - 1], pos[2*j]
                dx = x - jx
                dy = y - jy
                dist_sq = dx^2 + dy^2 + 1e-6
                
                force[idx_x] += k_p * dx / dist_sq
                force[idx_y] += k_p * dy / dist_sq
            end
        end
    end
    return force
end

function bounding_box_force!(force, pos, n_pts, bounds, k_wall)
    (xmin, xmax), (ymin, ymax) = bounds
    
    for i in 1:n_pts
        idx_x, idx_y = 2*i - 1, 2*i
        x, y = pos[idx_x], pos[idx_y]
        
        force[idx_x] += k_wall / (x - xmin)^2
        force[idx_x] -= k_wall / (xmax - x)^2
        
        force[idx_y] += k_wall / (y - ymin)^2
        force[idx_y] -= k_wall / (ymax - y)^2
    end
    return force
end

function torque_force!(force, pos, collinearity_relns, k_torque; max_cot=10.0)
    n_lines = length(collinearity_relns)
    
    for L1_idx in 1:n_lines
        for L2_idx in (L1_idx+1):n_lines
            L1 = collinearity_relns[L1_idx]
            L2 = collinearity_relns[L2_idx]
            
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
                        
                        # Tangential force on point i (removed division by norm2_A)
                        fti = -k_torque * cot_theta 
                        Fix = fti * (-Ay) 
                        Fiy = fti * (Ax)
                        force[idx_ix] += Fix
                        force[idx_iy] += Fiy
                        
                        # Equal and opposite tangential force on point j (removed division by norm2_B)
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