""" 
    AbstractObservables

Defines an interface to measure arbitrary observables during Monte Carlo simulations.
Observables can be provided by defining a new struct:
struct YourObservablesType <: AbstractObservables
    # Your fields here
end
and overloading Carlo.measure!(mc::AbstractWavefunctionMC, ctx::MCContext, obs::YourObservablesType).

Alternatively, you can also create a BufferedObservables instance, which supports
buffering intermediate results (provided via `buffer` and `buffer_functions`). E.g.
one can use this to calculate the pairwise distances between particles and then
use these distances to compute observables.
The observables themselves can be provided as name => function or as a struct <= Observables.
As a function, the signature should be
    f(bare_position, position, buffer) -> value
where bare_position is the raw position of the particles, position is the transformed
position (after applying the coordinate_transformer) and buffer is a tuple of all
buffers.

As a struct, you need to overload
    Carlo.measure!(mc::AbstractWavefunctionMC, ctx::MCContext, obs::YourObservableStruct, bare_position, position, buffer).

"""
abstract type AbstractObservables end
abstract type Observable end

struct NoObservables <: AbstractObservables end

@kwdef struct BufferedObservables{B <: Tuple, F <: Tuple, OF <: Tuple, OS <: Tuple} <: AbstractObservables
    buffer::B = NTuple{0, AbstractVector}()
    buffer_functions::F = NTuple{0, Function}()
    observable_functions::OF = NTuple{0, Pair{Symbol, Function}}()
    observable_structs::OS = NTuple{0, Observable}()
end

function distance!(buf::AbstractVector, pos, f::Function = Base.abs)
    N = length(pos)
    @boundscheck length(buf) >= div(N * (N - 1), 2)
    counter = 1
    @inbounds for i in 1:(N - 1)
        x = pos[i]
        for j in (i + 1):N
            buf[counter] = f(x - pos[j])
            counter += 1
        end
    end
    return buf
end

function Carlo.measure!(::AbstractWavefunctionMC, ::MCContext, ::NoObservables)
    return nothing
end

function Carlo.measure!(mc::MC, ctx::MCContext, obs::BufferedObservables) where {MC <: AbstractWavefunctionMC}
    bare_position = state(mc).bare_position
    pos = position(state(mc))
    buffer = obs.buffer
    for (b, f!) in zip(obs.buffer, obs.buffer_functions)
        f!(b, bare_position, pos)
    end
    for (name, f) in obs.observable_functions
        val = f(bare_position, pos, buffer)
        measure!(ctx, name, val)
    end
    for s in obs.observable_structs
        measure!(mc, ctx, s, bare_position, pos, buffer)
    end
    return nothing
end
