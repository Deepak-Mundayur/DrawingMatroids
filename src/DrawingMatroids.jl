module DrawingMatroids

import Random
using Random: rand, randn
import HomotopyContinuation
using HomotopyContinuation: System, Variable, expressions, variables
using Plots: Animation, annotate!, frame, gif, plot, plot!, savefig, scatter!, text
# import Oscar
using Oscar
using LinearAlgebra: I, det
using Printf: @sprintf

import Constrained_Optimization
const constrained_optimize = Constrained_Optimization.optimize


include(joinpath(@__DIR__, "..", "NumericalRealizationSpaces", "RealizationSpaces", "src", "RealizationSpaces.jl"))
include(joinpath(@__DIR__, "remi", "maximat.jl"))
import .Maximat

# using RealizationSpaces

# include("varieties.jl")
# include("vector_fields.jl")
# include("api.jl")

include("types.jl")
include("matroid_adapters.jl")
include("systems.jl")
include("initializers.jl")
include("validation.jl")
include("rendering.jl")
include("visualization.jl")
include("display.jl")


export visualization, matroid_visualize

export DrawingResult,
       DrawingCollectionResult,
       DrawingValidation,
       InitializationResult

export validation_summary,
       drawing_reduction

export maximal_degeneration_basis_transversal

# export matroid_collinearity_system, drawing_system,
#        collinearity_triples,
#        matroid_collinearity_vector_field,
#        matroid_collinearity_system_and_vector_field,

export render_drawing,
       save_static_drawing,
       save_drawing_animation,
       inspect

export constrained_optimize



end