abstract type AbstractAcceptanceAdapter end

struct NoAcceptanceAdapter <: AbstractAcceptanceAdapter end

@kwdef struct AcceptanceAdapter{T, F, G} <: AbstractAcceptanceAdapter
    acceptance_rate::T
    adapt_interval::Int

    # Determines, whether the adapter is only active during thermalization
    # or the entire run
    active_during_run::Bool = false

    # Defines the feedback function to adjust the acceptance rate
    # Should be > 1 for positive inputs and < 1 for negative inputs
    update_distribution::F = (d, adjustment) -> begin
        sigma_new = d.σ * Base.exp(-adjustment)
        return typeof(d)(d.μ, sigma_new)
    end
    update::G = Base.exp
end
