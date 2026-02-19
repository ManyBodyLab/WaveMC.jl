""" 
    WavefunctionMC 

Main type to run wavefunction-based Monte Carlo simulations.
"""
struct WavefunctionMC{N, S<:State, D <: Distribution, F, X <: AbstractWavefunction, obs <: AbstractObservables, A <: AbstractAcceptanceAdapter, B, U<:AbstractCoordinateUpdate} <: AbstractWavefunctionMC{N}
    state::S
    distribution::Ref{D}
    coordinate_projector::F ## Makes sure the positions are in the right interval
    wavefunction::X
    observables::obs   ## They can e.g. have a buffer inside for optimization!
    acceptance_adapter::A
    dynamic_positions::B
    coordinate_update::U
    function WavefunctionMC{N}(
            state::S,
            distribution::D,
            coordinate_projector::F,
            wavefunction::X,
            observables::obs;
            adaptive::Bool = true,
            adapt_interval::Int = 1_000,
            coordinate_update::U = CoordinateUpdate(),
            target_accept::Float64 = optimal_acceptance_rate(coordinate_update), 
            dynamic_positions::B = (N, 1:N),
            adaptor::A = adaptive ? AcceptanceAdapter(; acceptance_rate = target_accept, adapt_interval = adapt_interval) : NoAcceptanceAdapter()
        ) where {N, S, D <: Distribution, F, X <: AbstractWavefunction, obs <: AbstractObservables, B, A, U}
        return new{N, S, D, F, X, obs, A, B, U}(state, Ref(distribution), coordinate_projector, wavefunction, observables, adaptor, dynamic_positions, coordinate_update)
    end
end

@inline function state(mc::WavefunctionMC{N, S})::S where {N, S} 
    return mc.state
end
@inline function distribution(mc::WavefunctionMC{N, S, D})::D where {N, S, D <: Distribution} 
    return mc.distribution[]
end
@inline coordinate_projector(mc::WavefunctionMC) = mc.coordinate_projector
@inline wavefunction(mc::WavefunctionMC) = mc.wavefunction
@inline evaluables(mc::WavefunctionMC) = mc.evaluables
@inline acceptance_rate_adapter(mc::WavefunctionMC) = mc.acceptance_adapter
@inline dynamic_positions(mc::WavefunctionMC) = mc.dynamic_positions
@inline coordinate_update(mc::WavefunctionMC) = mc.coordinate_update
@inline observables(mc::WavefunctionMC) = mc.observables

function WavefunctionMC(params::AbstractDict)
    wavefunction = params[:wavefunction]

    coordinate_proj = haskey(params, :coordinate_projector) ? params[:coordinate_projector] : coordinate_projector(wavefunction)
    coordinate_transf = haskey(params, :coordinate_transformer) ? params[:coordinate_transformer] : coordinate_transformer(wavefunction)
    S = inputtype(wavefunction)
    N = inputlength(wavefunction)
    sigma_dist = get(params, :sigma_distribution, 0.3)
    distribution = get(params, :distribution, S <: Complex ? ComplexNormal(0, sigma_dist) : Normal(0, sigma_dist))
    position = coordinate_proj.(100 * rand(distribution, N))

    observables = get(params, :observables, NoObservables())
    dynamic_pos = get(params, :dynamic_positions, (N, 1:N))
    return WavefunctionMC{N}(State(position, 0.0, coordinate_transf), distribution, coordinate_proj, wavefunction, observables; adapt_interval = div(get(params, :thermalization, 10_000), 10), dynamic_positions = dynamic_pos)
end

function Carlo.write_checkpoint(mc::WavefunctionMC, out::HDF5.Group)
    out["pos_buf"] = mc.state.position.buffer
    out["logdensity"] = mc.state.logdensity
    out["num_accepts"] = mc.state.num_accepts
    out["num_total"] = mc.state.num_total
    return nothing
end

function Carlo.read_checkpoint!(mc::WavefunctionMC, in::HDF5.Group)
    mc_state = state(mc)
    mc_state.position.buffer .= read(in, "pos_buf")
    mc_state.logdensity = read(in, "logdensity")
    mc_state.num_accepts = read(in, "num_accepts")
    mc_state.num_total = read(in, "num_total")
    return nothing
end
