# function matroid_visualize(M)

#     return nothing
# end



"""
    run_and_animate_collinear_system(n_pts, bounds, collinear_sets; kwargs...)

Creates the matroid geometry, defines the force field, finds an equilibrium, and outputs a GIF.
"""
function run_and_animate_collinear_system(n_pts::Int, collinear_sets; 
                                          bounds = ((-5.0,5.0),(-5.0,5.0)),
                                          dt=0.01, max_steps=300, 
                                          k_p=10.0, k_wall=50.0,
                                          guess = nothing,
                                          filename="example.gif", fps=15,
                                          show_guess = false,
                                          framestyle = :none,
                                          hard_non_collinearity_relns = false,
                                          )
    
    (x_min, x_max), (y_min, y_max) = bounds

    
    F_sys = matroid_collinearity_system(n_pts, collinear_sets; hard_non_collinearity_relns=hard_non_collinearity_relns)
    n_vars = length(variables(F_sys))
    
    V = make_bounded_repelling_force(n_pts, bounds; k_p=k_p, k_wall=k_wall)
    
    if guess == nothing
        guess = Float64[]
        
        radius = min(x_max - x_min, y_max - y_min) * 0.35
        center_x = (x_max + x_min) / 2.0
        center_y = (y_max + y_min) / 2.0
        
        for i in 1:n_pts
            # Distribute points evenly around a circle
            angle = 2 * pi * i / n_pts
            push!(guess, center_x + radius * cos(angle))
            push!(guess, center_y + radius * sin(angle))
        end
        
        # Fill remaining slack variables if using the hard-boundary matroid system
        num_slacks = n_vars - length(guess)
        if num_slacks > 0
            append!(guess, fill(0.5, num_slacks))
        end
    end    

    guess_xs = [guess[2*i - 1] for i in 1:n_pts]
    guess_ys = [guess[2*i] for i in 1:n_pts]
    
    final_point, history = constrained_optimize(F_sys, V; guess=guess, dt=dt, max_steps=max_steps, corrector = moore_penrose_corrector)
    anim = Animation()
    
    for (step_idx, state) in enumerate(history)
        xs = [state[2*i - 1] for i in 1:n_pts]
        ys = [state[2*i] for i in 1:n_pts]
        
        plt = plot(
            xlim=(x_min - 1.0, x_max + 1.0), 
            ylim=(y_min - 1.0, y_max + 1.0),
            aspect_ratio=:equal,
            legend=false,
            title="DrawingMatroids.jl",
            framestyle = framestyle
        )
        
        plot!(plt, [x_min, x_max, x_max, x_min, x_min], 
                   [y_min, y_min, y_max, y_max, y_min], 
                   color=:red, linewidth=2, linestyle=:dash)
        
        for (p1, p2, p3) in collinear_sets
            dx = xs[p2] - xs[p1]
            dy = ys[p2] - ys[p1]
            
            # Safety check: If p1 and p2 are completely crushed together 
            # (like in step 1), use p3 to calculate the line instead
            if abs(dx) < 1e-5 && abs(dy) < 1e-5
                dx = xs[p3] - xs[p1]
                dy = ys[p3] - ys[p1]
            end
            
            # If all 3 points are crushed into the exact same pixel, skip the line for this frame
            if abs(dx) < 1e-5 && abs(dy) < 1e-5
                continue
            end

            if abs(dx) > 1e-7
                slope = dy / dx
                intercept = ys[p1] - slope * xs[p1]
                line_xs = [x_min - 2.0, x_max + 2.0]
                line_ys = slope .* line_xs .+ intercept
                plot!(plt, line_xs, line_ys, color=:gray, alpha=0.4, linewidth=1.5)
            else
                # Handle perfectly vertical lines
                plot!(plt, [xs[p1], xs[p1]], [y_min - 2.0, y_max + 2.0], color=:gray, alpha=0.4, linewidth=1.5)
            end
        end

        # Scatter the initial guess in the background
        if show_guess
            scatter!(plt, guess_xs, guess_ys, color=:lightgray, markersize=8, alpha=0.5)
        end 

        # Scatter the active beads
        scatter!(plt, xs, ys, color=:blue, markersize=8)
        
        frame(anim, plt)
    end
    
    gif(anim, filename, fps=fps)
    println("Animation saved to: $filename")
    
    return final_point, history
end


function run_and_animate_collinear_system(n_pts::Int, bounds, collinear_sets; kwargs...)
    return run_and_animate_collinear_system(n_pts, collinear_sets; bounds = bounds, kwargs...)
end
