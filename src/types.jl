struct NotComplexRealizableError <: Exception
    message::String
end

Base.showerror(io::IO, err::NotComplexRealizableError) = print(io, err.message)

struct NoRealDrawingPointError <: Exception
    message::String
end

Base.showerror(io::IO, err::NoRealDrawingPointError) = print(io, err.message)

struct RelabeledMatroidData
    original
    matroid
    labels::Vector{Any}
    label_to_index::Dict{Any,Int}
end

struct DrawingReduction
    relabeled::RelabeledMatroidData
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
