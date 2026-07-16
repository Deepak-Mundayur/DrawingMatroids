function matroid_collinearity_system(n_points::Int, collinear_sets::Vector{NTuple{3, Int}}; hard_non_collinearity_relns=false)
    # Variables for n points (x, y)
    vars_xy = Variable[]
    for i in 1:n_points
        push!(vars_xy, Variable("x$i"))
        push!(vars_xy, Variable("y$i"))
    end
    
    eqs = Expression[]
    
    # Sorting
    collinear_lookup = Set(Tuple(sort([c...])) for c in collinear_sets)
    
    # Collinearity Equations
    for (i, j, k) in collinear_sets
        xi, yi = vars_xy[2*i - 1], vars_xy[2*i]
        xj, yj = vars_xy[2*j - 1], vars_xy[2*j]
        xk, yk = vars_xy[2*k - 1], vars_xy[2*k]
        
        collinear_eq = (xj - xi)*(yk - yi) - (yj - yi)*(xk - xi)
        push!(eqs, collinear_eq)
    end
    
    # The non-collinearity conditions
    t_vars = Variable[]
    
    if hard_non_collinearity_relns
        # Iterate over all unique combinations of 3 points
        for i in 1:(n_points - 2)
            for j in (i + 1):(n_points - 1)
                for k in (j + 1):n_points
                    
                    # If this triplet is NOT meant to be collinear, force the determinant away from 0
                    if !((i, j, k) in collinear_lookup)
                        # Create a unique slack variable for this specific non-collinear triplet
                        t_var = Variable("t_$(i)_$(j)_$(k)")
                        push!(t_vars, t_var)
                        
                        xi, yi = vars_xy[2*i - 1], vars_xy[2*i]
                        xj, yj = vars_xy[2*j - 1], vars_xy[2*j]
                        xk, yk = vars_xy[2*k - 1], vars_xy[2*k]
                        
                        det_ijk = (xj - xi)*(yk - yi) - (yj - yi)*(xk - xi)
                        push!(eqs, t_var * det_ijk - 1.0)
                    end
                end
            end
        end
        return System(eqs, variables=vcat(vars_xy, t_vars))
    end
    
    return System(eqs, variables=vars_xy)
end

