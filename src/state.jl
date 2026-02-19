mutable struct State{S, B <: Buffer, F, C <: Integer} 
    bare_position::Vector{S}
    position::B
    logdensity::F
    num_accepts::C
    num_total::C
end

position(state::State) = state.position
function Base.getindex(state::State, i)
    return state.position[i]
end
function Base.setindex!(state::State, v, i)
    state.bare_position[i] = v
    return state.position[i] = v
end

function coordinate_transformer(state::State)
    return coordinate_transformer(state.position)
end
function State(physical_position::AbstractVector{T}, logdensity::Float64, coordinate_transformer::G = Base.identity) where {T, G}
    pos = coordinate_transformer.(physical_position)
    #pos = MVector{length(physical_position), eltype(pos)}(pos)
    return State(physical_position, Buffer(copy(pos); map = coordinate_transformer), logdensity, 0, 0)
end

function acceptance_rate(state::State)
    state.num_total == 0 && return 1.0
    return state.num_accepts / state.num_total
end

function reset_acceptance_counters!(state::State)
    state.num_accepts = 0
    state.num_total = 0
    return nothing
end
