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
    (xmin, xmax), (ymin, ymax) = bounds

    collinearity_relns = collect.(collinearity_tuples)

    function V(pos)
        force = zeros(length(pos))
        
        # 1. Point-to-Point and Bounding Box Forces
        for i in 1:n_pts
            idx_x = 2*i - 1
            idx_y = 2*i
            x, y = pos[idx_x], pos[idx_y]
            
            # Point-to-Point Repulsion 
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
            
            # The bounding box forces (Inverse square barrier)
            force[idx_x] += k_wall / (x - xmin)^2
            force[idx_x] -= k_wall / (xmax - x)^2
            
            force[idx_y] += k_wall / (y - ymin)^2
            force[idx_y] -= k_wall / (ymax - y)^2
        end
        
        # 2. Line-to-Line Torque Repulsion
        n_lines = length(collinearity_relns)
        for L1_idx in 1:n_lines
            for L2_idx in (L1_idx+1):n_lines
                L1 = collinearity_relns[L1_idx]
                L2 = collinearity_relns[L2_idx]
                
                # Find the intersection node (center of torque)
                shared_pts = intersect(L1, L2)
                if length(shared_pts) == 1
                    c = shared_pts[1]
                    cx, cy = pos[2*c - 1], pos[2*c]
                    
                    # Apply torque between all pairs of points across the two lines
                    for i in L1
                        if i == c continue end
                        idx_ix, idx_iy = 2*i - 1, 2*i
                        Ax, Ay = pos[idx_ix] - cx, pos[idx_iy] - cy
                        norm2_A = Ax^2 + Ay^2
                        if norm2_A < 1e-6 continue end # Safety check
                        
                        for j in L2
                            if j == c continue end
                            idx_jx, idx_jy = 2*j - 1, 2*j
                            Bx, By = pos[idx_jx] - cx, pos[idx_jy] - cy
                            norm2_B = Bx^2 + By^2
                            if norm2_B < 1e-6 continue end # Safety check
                            
                            # Dot product (proportional to cos θ) and Cross product (proportional to sin θ)
                            dp = Ax * Bx + Ay * By
                            cp = Ax * By - Ay * Bx
                            
                            # Prevent division by zero if lines perfectly overlap
                            cp_safe = cp >= 0.0 ? max(cp, 1e-6) : min(cp, -1e-6)
                            
                            # The cotangent (dp/cp) dictates the magnitude and direction of the angular force.
                            # Dividing by norm2 creates a physical force that acts like a true torque (F = τ/r)
                            
                            # Tangential force on point i
                            fti = -k_torque * (dp / cp_safe) / norm2_A
                            force[idx_ix] += fti * (-Ay) # -Ay, Ax is the tangential CCW direction
                            force[idx_iy] += fti * (Ax)
                            
                            # Equal and opposite tangential force on point j
                            ftj = k_torque * (dp / cp_safe) / norm2_B
                            force[idx_jx] += ftj * (-By)
                            force[idx_jy] += ftj * (Bx)
                        end
                    end
                end
            end
        end
        
        # If the system has a saturation slack variable at the end, we assume that the its component of the vector field in that direction is 0 
        return force
    end
    
    return V
end