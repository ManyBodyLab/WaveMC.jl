abstract type AbstractCoordinateUpdate end

struct GlobalUpdate <: AbstractCoordinateUpdate end 
struct CoordinateUpdate <: AbstractCoordinateUpdate end

function optimal_acceptance_rate(x::T) where {T<:AbstractCoordinateUpdate} 
    return optimal_acceptance_rate(T)
end
optimal_acceptance_rate(::Type{AbstractCoordinateUpdate}) = 0.234
optimal_acceptance_rate(::Type{CoordinateUpdate}) = 0.44
