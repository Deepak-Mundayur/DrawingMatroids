module DrawingMatroids

import Random
using Random: rand, randn
using LinearAlgebra: I, det, svd

import HomotopyContinuation
using HomotopyContinuation: System, Expression, @var
using Plots: Animation, annotate!, frame, gif, plot, plot!, savefig, scatter!, text
using Oscar

import Constrained_Optimization: optimize as constrained_optimize, moore_penrose_corrector, project_onto_manifold, make_bounded_repelling_force
using RealizationSpaces
include(joinpath(@__DIR__, "maximat", "maximat.jl"))
import .Maximat

include("matroid_adapters.jl")
include("systems.jl")
include("initializers.jl")
include("validation.jl")
include("rendering.jl")
include("vector_fields.jl")
include("visualization.jl")


export draw_matroid, draw_component

end