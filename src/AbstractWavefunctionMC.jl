""" 
    AbstractWavefunctionMC{N}

Abstract type for wavefunction-based Monte Carlo methods with N coordinates.

Functions to be implemented for concrete subtypes:
- state(mc::AbstractWavefunctionMC): Returns the state of the Monte Carlo simulation.
- distribution(mc::AbstractWavefunctionMC): Returns the proposal distribution used for updates.
- wavefunction(mc::AbstractWavefunctionMC): Returns the wavefunction used in the simulation.
- observables(mc::AbstractWavefunctionMC): Returns the observables to be measured.
- dynamic_positions(mc::AbstractWavefunctionMC{N}): Returns a tuple (N_active, dynamic_pos) indicating the number of active positions and their indices.
- acceptance_rate_adapter(mc::AbstractWavefunctionMC): Returns the acceptance rate adapter used for adapting the proposal distribution. (defaults to NoAcceptanceAdapter)
- change_distribution!(mc::AbstractWavefunctionMC, new_dist::Distribution): Changes the proposal distribution to new_dist. (Only necessary for acceptance_rate_adapter)
"""
abstract type AbstractWavefunctionMC{N} <: AbstractMC end

evaluables(::AbstractWavefunctionMC) = nothing
dynamic_positions(::AbstractWavefunctionMC{N}) where {N} = (N, 1:N)
function acceptance_rate_adapter(mc::AbstractWavefunctionMC)
    return NoAcceptanceAdapter()
end
function change_distribution!(mc::AbstractWavefunctionMC, new_dist::Distribution)
    mc.distribution[] = new_dist
    return nothing
end

# --------- Utility functions & Carlo.jl ---------
function coordinate_transformer(mc::AbstractWavefunctionMC)
    return coordinate_transformer(state(mc))
end
acceptance_rate(mc::AbstractWavefunctionMC) = acceptance_rate(state(mc))

function Carlo.init!(mc::AbstractWavefunctionMC{N}, ctx::MCContext, params::AbstractDict) where {N}
    mc_state = state(mc)
    N_active, dynamic_pos = dynamic_positions(mc)

    mc_state[:] = coordinate_projector(mc).(rand(ctx.rng, distribution(mc), N))

    mc_state.logdensity = logdensity(wavefunction(mc), position(mc_state).buffer)
    mc_state.num_accepts = 0
    mc_state.num_total = 0
    return nothing
end

function Carlo.sweep!(mc::MC, ctx::MCContext) where {MC <: AbstractWavefunctionMC{N}} where {N}
    adapt_acceptance_rate!(mc, ctx)
    return Carlo.sweep!(mc, ctx, coordinate_update(mc))
end
function Carlo.sweep!(mc::MC, ctx::MCContext, ::CoordinateUpdate) where {MC <: AbstractWavefunctionMC{N}} where {N}
    d = distribution(mc)
    mc_state = state(mc)
    bare_position, buffer = mc_state.bare_position, mc_state.position

    coordinate_trafo = coordinate_transformer(mc)
    coordinate_proj = coordinate_projector(mc)
    wavefunc = wavefunction(mc)
    num_updates, dynamic_pos = dynamic_positions(mc)
    _update_coordinates(ctx, dynamic_pos, bare_position, d, coordinate_trafo, coordinate_proj, wavefunc, buffer, mc_state)
    mc_state.num_total += num_updates
    return nothing
end

function _update_coordinates(ctx::MCContext, dynamic_pos, bare_position, d, coordinate_trafo, coordinate_proj, wavefunc::T, buffer::X, mc_state) where {T<:AbstractWavefunction, X<:Buffer}
    @inbounds for dim in dynamic_pos
        x_new = coordinate_proj(bare_position[dim] + rand(ctx.rng, d))
        x_new_transformed = coordinate_trafo(x_new)

        logα = delta_logdensity(wavefunc, x_new_transformed, buffer.buffer, dim)
        if -randexp(ctx.rng) < logα
            logdens_new = mc_state.logdensity + logα
            bare_position[dim] = x_new
            buffer.buffer[dim] = x_new_transformed
            mc_state.logdensity = logdens_new
            mc_state.num_accepts += 1
        end
    end
    return nothing
end

function Carlo.sweep!(mc::MC, ctx::MCContext, ::GlobalUpdate) where {MC <: AbstractWavefunctionMC{N}} where {N}
    d = distribution(mc)
    mc_state = state(mc)
    bare_position, buffer = mc_state.bare_position, mc_state.position

    coordinate_trafo = coordinate_transformer(mc)
    wavefunc = wavefunction(mc)
    _, dynamic_pos = dynamic_positions(mc)

    pos_new = mc_state.bare_position .+ rand(ctx.rng, d, N)
    pos_new_transformed = coordinate_trafo.(pos_new)
    logα = delta_logdensity(wavefunc, pos_new_transformed, buffer.buffer)
    if -randexp(ctx.rng) < logα
        logdens_new = mc_state.logdensity + logα
        for dim in dynamic_pos
            bare_position[dim] = pos_new[dim]
            buffer.buffer_previous[dim] = buffer.buffer[dim]
            buffer.buffer[dim] = pos_new_transformed[dim]
        end
        mc_state.logdensity = logdens_new
        mc_state.num_accepts += 1
    end
    mc_state.num_total += 1
    return nothing
end

@inline function Carlo.measure!(mc::AbstractWavefunctionMC, ctx::MCContext)
    measure!(ctx, :acceptance_rate, acceptance_rate(mc))
    measure!(mc, ctx, observables(mc))
    return nothing
end

# TODO: Support evaluables!
function Carlo.register_evaluables(
        ::Type{<:AbstractWavefunctionMC},
        eval::AbstractEvaluator,
        params::AbstractDict,
    )
    return nothing
end
