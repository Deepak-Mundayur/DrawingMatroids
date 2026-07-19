module DrawingMatroids

using Random
using HomotopyContinuation
using Plots
using Oscar
using LinearAlgebra
using Printf

import Constrained_Optimization
const constrained_optimize = Constrained_Optimization.optimize
const moore_penrose_corrector = Constrained_Optimization.moore_penrose_corrector

include(joinpath(@__DIR__, "..", "NumericalRealizationSpaces", "RealizationSpaces", "src", "RealizationSpaces.jl"))
include(joinpath(@__DIR__, "remi", "maximat.jl"))
import .Maximat

# using RealizationSpaces

# include("varieties.jl")
# include("vector_fields.jl")

include("types.jl")
include("matroid_adapters.jl")
include("systems.jl")
include("initializers.jl")
include("validation.jl")
include("rendering.jl")
include("visualization.jl")
include("api.jl")


export RealizationSpaces, Maximat

export DrawingResult,
       DrawingCollectionResult,
       DrawingValidation,
       InitializationResult,
       NotComplexRealizableError,
       NoRealDrawingPointError

export visualization,
       draw,
       matroid_visualize,
       validation_summary

export drawing_reduction,
       relabel_matroid,
       realization_matrix_to_coordinates,
       coordinates_to_matrix,
       cyclic_flats_by_rank,
       maximal_degeneration_basis_transversal,
       rank_two_flats

export matroid_collinearity_system,
       drawing_system,
       collinearity_triples,
       circle_guess,
       append_rabinowitsch_coordinates

export render_drawing,
       save_static_drawing,
       save_drawing_animation,
       run_and_animate_collinear_system

export constrained_optimize

# export matroid_collinearity_system, make_bounded_repelling_force #matroid_visualize
# export run_and_animate_collinear_system

end