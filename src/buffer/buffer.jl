"""   
    Buffer{T, A<:AbstractVector{T}, G} 
is a wrapper around a AbstractVector{T} that applies a mapping G to values before storing them in the buffer.

Fields:
- buffer::A: The underlying dense vector buffer.
- map::G: A function that transforms values before they are stored in the buffer. Defaults to the identity function.

This structure allows for efficient storage and retrieval of values with an optional transformation applied during storage.
"""
@kwdef struct Buffer{T, A <: AbstractVector{T}, G} <: AbstractVector{T}
    buffer::A
    map::G = Base.identity
end
Base.parent(b::Buffer) = b.buffer
Base.size(b::Buffer) = size(parent(b))
Base.length(b::Buffer) = length(parent(b))
Base.@propagate_inbounds Base.getindex(b::Buffer, i) = getindex(parent(b), i)
Base.firstindex(b::Buffer) = firstindex(parent(b))
Base.lastindex(b::Buffer) = lastindex(parent(b))
Base.@propagate_inbounds Base.setindex!(b::Buffer, v, i) = (parent(b)[i] = b.map(v))

function Buffer(def::T, size; map::G = Base.identity) where {T, G}
    return Buffer{T, AbstractVector{T}, G}(fill(def, size), map)
end
function Buffer(vec::AbstractVector{T}; map::G = Base.identity) where {T, G}
    return Buffer{T, AbstractVector{T}, G}(vec, map)
end

function coordinate_transformer(b::Buffer)
    return b.map
end

# """
#     Buffer{T, A<:AbstractVector{T}, G} 
# extends Buffer by storing the previous values before they are updated.

# Fields:
# - buffer::Buffer{T, A, G}: The underlying buffer that stores current values.
# - buffer_previous::A: A vector that stores the previous values before updates.

# This structure is useful for scenarios where tracking changes to buffer values is necessary.
# """
# @kwdef struct Buffer{T, A <: AbstractVector{T}, G}
#     buffer::Buffer{T, A, G}
#     buffer_previous::A
# end

# Base.parent(b::Buffer) = b.buffer
# Base.size(b::Buffer) = size(parent(b))
# Base.length(b::Buffer) = length(parent(b))
# Base.@propagate_inbounds Base.getindex(b::Buffer, i) = getindex(parent(b), i)
# Base.firstindex(b::Buffer) = firstindex(parent(b))
# Base.lastindex(b::Buffer) = lastindex(parent(b))
# Base.@propagate_inbounds function Base.setindex!(b::Buffer, v, i)
#     b.buffer_previous[i] = parent(b)[i]
#     return parent(b)[i] = v
# end

# function Buffer(def::T, size; map::G = Base.identity) where {T, G}
#     buf = Buffer{T, AbstractVector{T}, G}(fill(def, size), map)
#     buf_prev = copy(parent(buf))
#     return Buffer{T, AbstractVector{T}, G}(buf, buf_prev)
# end
# function Buffer(vec::AbstractVector{T}; map::G = Base.identity) where {T, G}
#     buf = Buffer{T, AbstractVector{T}, G}(vec, map)
#     buf_prev = copy(parent(buf))
#     return Buffer{T, AbstractVector{T}, G}(buf, buf_prev)
# end

# function coordinate_transformer(b::Buffer)
#     return coordinate_transformer(b.buffer)
# end
