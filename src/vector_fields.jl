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