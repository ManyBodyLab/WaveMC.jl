abstract type AbstractBuffer end

function fill_buffer!(buf::AbstractBuffer, args...)
    return buf
end

struct NoBuffer <: AbstractBuffer 
    NoBuffer(args...; kwargs...) = new() 
end

