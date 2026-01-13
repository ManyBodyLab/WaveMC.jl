
function adapt_acceptance_rate!(psi::AbstractWavefunctionMC, ctx::MCContext)
    return adapt_acceptance_rate!(psi, ctx, acceptance_rate_adapter(psi))
end
function adapt_acceptance_rate!(::AbstractWavefunctionMC, ::MCContext, ::NoAcceptanceAdapter)
    mc_state = state(mc)
    return reset_acceptance_counters!(mc_state)
end
function adapt_acceptance_rate!(mc::MC, ctx::MCContext, adapter::GaussianAcceptanceAdapter; verbosity::Int=1) where {MC<:AbstractWavefunctionMC}
    mc_state = state(mc)
    if !is_thermalized(ctx) && mod1(ctx.sweeps, adapter.adapt_interval) == adapter.adapt_interval
        acc_rate = acceptance_rate(mc)

        if verbosity > 0
            str = sprint("Adapt acceptance rate: Current acceptance rate = %.4f, target = %.4f", acc_rate, adapter.acceptance_rate)
            @info str 
        end

        adjustment = if acc_rate < adapter.acceptance_rate 
            adapter.update(acc_rate - adapter.acceptance_rate)
        else
            adapter.update(acc_rate - adapter.acceptance_rate)
        end

        d = distribution(mc)
        sigma_new = d.σ * adjustment
        mc.distribution[] = typeof(d)(d.μ, sigma_new)
    end

    return reset_acceptance_counters!(mc_state)
end
