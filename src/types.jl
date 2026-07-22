struct DrawingReduction
    matroid                              # relabeled matroid on ground set 1:n
    simple_matroid
    parallel_classes::Vector{Vector{Int}}
    representatives::Vector{Int}
    grouped_labels::Vector{Vector{Any}}
    loops::Vector{Int}
    loop_labels::Vector{Any}
end

struct InitializationResult
    strategy::Symbol
    state::Vector{Float64}
    system::System
    separating_bases::Vector{Vector{Int}}
    maximal_degenerations::Vector{Any}
    component_index::Union{Nothing,Int}
    diagnostics::Dict{Symbol,Any}
end

struct DrawingValidation
    valid::Bool
    matrix::Matrix{Float64}
    maximum_nonbasis_residual::Float64
    minimum_basis_absolute_determinant::Float64
    failed_nonbases::Vector{Vector{Int}}
    failed_bases::Vector{Vector{Int}}
    collisions::Vector{Tuple{Int,Int}}
    component_member::Union{Nothing,Bool}
end

struct DrawingResult
    source
    target_matroid
    reduction::DrawingReduction
    initialization::InitializationResult
    final_state::Vector{Float64}
    history::Vector{Vector{Float64}}
    validation::DrawingValidation
    plot
    image_filename::Union{Nothing,String}
    animation_filename::Union{Nothing,String}
end

struct DrawingCollectionResult
    results::Vector{DrawingResult}
    plot
    filename::Union{Nothing,String}
end

Base.@kwdef struct DrawingOptions
    # Constrained optimization
    dt::Float64 = 0.01
    max_steps::Int = 300
    k_p::Float64 = 10.0
    k_wall::Float64 = 50.0

    # Initial projection
    projection_tol::Float64 = 1e-10
    projection_max_iters::Int = 100

    # Numerical validation
    validation_tol::Float64 = 1e-7
    collision_tol::Float64 = 1e-6
    real_tol::Float64 = 1e-7

    # Rémi
    remi_preprocess::Bool = true
    remi_verbosity::Int = 0
    remi_threads::Int = Base.Threads.nthreads()

    # Rendering
    fps::Int = 15
    title::String = "DrawingMatroids.jl"
    framestyle::Symbol = :box
    show_boundary::Bool = true
    show_labels::Bool = true
end