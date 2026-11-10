module DrawingMatroids

import Random
using Random: rand, randn
using LinearAlgebra: I, det, svd

import HomotopyContinuation
import JuMP, Ipopt
using HomotopyContinuation: System, Expression, @var
using Plots: Animation, Plot, Shape, annotate!, frame, gif, palette, plot, plot!, savefig, scatter!, text
using Oscar

import Constrained_Optimization: optimize as constrained_optimize, moore_penrose_corrector, project_onto_manifold, make_bounded_repelling_force
using RealizationSpaces
include(joinpath(@__DIR__, "maximat", "maximat.jl"))
import .Maximat

include("coordinates.jl")
include("systems.jl")
include("initializers.jl")
include("validation.jl")
include("rendering.jl")
include("readability.jl")
include("forces.jl")
include("visualization.jl")
include("poset.jl")


export draw_matroid, draw_poset

end