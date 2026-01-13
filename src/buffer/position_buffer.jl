
@kwdef struct PositionBuffer{T, A<:AbstractVector{T}, G}
    position_buffer::A
    coordinate_transformer::G=Base.identity
end

@kwdef struct Position2Buffer{T, A<:AbstractVector{T}, G}
    position_buffer_prev::A
    position_buffer::A = copy(position_buffer_prev)
    coordinate_transformer::G = Base.identity
end

function PositionBuffer(def::T, size, coordinate_transformer::G=Base.identity) where {T; G}
    return PositionBuffer(fill(def, size), coordinate_transformer)
end
function Position2Buffer(N::Int, T::Type=Float64; coordinate_transformer::G=Base.identity) where {G}
    pos_buf_prev = zeros(T, N)
    return Position2Buffer{T,G}(;pos_buf_prev = pos_buf_prev, coordinate_transformer=coordinate_transformer)
end

function fill_buffer!(f::PositionBuffer{T,G}, x::AbstractVector{T2}) where {T,G,T2<:Number}
    @inbounds @simd for i in eachindex(x)
        fill_buffer!(f, x[i], i)
    end
    return nothing
end
function fill_buffer!(f::PositionBuffer{T,G}, x::AbstractVector{T2}, ind::Int) where {T,G,T2<:Number}
    return fill_buffer!(f, x[ind], ind)
end
function fill_buffer!(f::PositionBuffer{T,G}, x::T2, ind::Int) where {T,G,T2<:Number}
    @inbounds f.pos_buf[ind] = f.coordinate_transformer(x)
    return nothing
end

function fill_buffer!(f::Position2Buffer{T,G}, x::AbstractVector{T2}) where {T,G,T2<:Number}
    @inbounds @simd for i in eachindex(x)
        fill_buffer!(f, x[i], i)
    end
    return nothing
end
@inline function fill_buffer!(f::Position2Buffer{T,G}, x::AbstractVector{T2}, ind::Int) where {T,G,T2<:Number}
    return fill_buffer!(f, x[ind], ind)
end
@inline function fill_buffer!(f::Position2Buffer{T,G}, x::T2, ind::Int) where {T,G,T2<:Number}
    @inbounds f.pos_buf_prev[ind] = f.pos_buf[ind]
    @inbounds f.pos_buf[ind] = f.coordinate_transformer(x)
    return nothing
end
@inline function fill_buffer_no_transform!(f::Position2Buffer{T,G}, x::T2, ind::Int) where {T,G,T2<:Number}
    @inbounds f.pos_buf_prev[ind] = f.pos_buf[ind]
    @inbounds f.pos_buf[ind] = x
    return nothing
end
function fill_buffer!(f::Position2Buffer{T,G}, xnew::AbstractVector{T2}, xold::AbstractVector{T2}) where {T,G,T2<:Number}
    @inbounds @simd for i in eachindex(xnew)
        f.pos_buf[i] = f.coordinate_transformer(xnew[i])
        f.pos_buf_prev[i] = f.coordinate_transformer(xold[i])
    end
    return nothing
end
