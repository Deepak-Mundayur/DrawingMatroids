module DrawingMatroids

using HomotopyContinuation
using Plots
using LinearAlgebra

#import Constrained_Optimization: optimize

include("varieties.jl")
include("vector_fields.jl")
include("visualization.jl")

export matroid_visualize

end