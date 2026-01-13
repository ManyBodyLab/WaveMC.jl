abstract type AbstractAcceptanceAdapter end 

struct NoAcceptanceAdapter <: AbstractAcceptanceAdapter end

@kwdef struct GaussianAcceptanceAdapter{T, F} <: AbstractAcceptanceAdapter
    acceptance_rate::T
    adapt_interval::Int

    # Determines, whether the adapter is only active during thermalization
    # or the entire run
    active_during_run::Bool = true

    # Defines the feedback function to adjust the acceptance rate
    # Should be > 1 for positive inputs and < 1 for negative inputs
    update::F = Base.exp    
end

