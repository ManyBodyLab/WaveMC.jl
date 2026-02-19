function adapt_acceptance_rate!(psi::AbstractWavefunctionMC, ctx::MCContext)
    return adapt_acceptance_rate!(psi, ctx, acceptance_rate_adapter(psi))
end
function adapt_acceptance_rate!(mc::AbstractWavefunctionMC, ::MCContext, ::NoAcceptanceAdapter)
    mc_state = state(mc)
    return reset_acceptance_counters!(mc_state)
end
function adapt_acceptance_rate!(mc::MC, ctx::MCContext, adapter::AcceptanceAdapter; verbosity::Int = 1) where {MC <: AbstractWavefunctionMC}
    mc_state = state(mc)
    if !is_thermalized(ctx) && mod1(ctx.sweeps, adapter.adapt_interval) == adapter.adapt_interval
        acc_rate = acceptance_rate(mc)
        if verbosity > 0
            str = @sprintf("Adapt acceptance rate: Current acceptance rate = %.6f, target = %.6f", acc_rate, adapter.acceptance_rate)
            @info str
        end

        change_distribution!(mc, adapter.update_distribution(distribution(mc), adapter.acceptance_rate - acc_rate))
        reset_acceptance_counters!(mc_state)
    end

    return nothing
end
