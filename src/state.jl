mutable struct State{S, X<:AbstractBuffer, F, C}
    bare_position::Vector{S}
    position::X
    logdensity::F
    num_accepts::C
    num_total::C
end
