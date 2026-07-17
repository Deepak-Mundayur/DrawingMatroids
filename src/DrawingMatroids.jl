module DrawingMatroids

using HomotopyContinuation
using Plots
using LinearAlgebra

import Constrained_Optimization: optimize as constrained_optimize, moore_penrose_corrector

include("varieties.jl")
include("vector_fields.jl")
include("visualization.jl")

export constrained_optimize
export matroid_collinearity_system, make_bounded_repelling_force #matroid_visualize
export run_and_animate_collinear_system

end