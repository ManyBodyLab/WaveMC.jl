struct Data{T, L<:Real, A<:AbstractVector{T}}
    data::A
    weight::L 
end 
const DataSet{T, L, A} = Vector{Data{T,L, A}} where {A<:AbstractVector{T}} where {T, L<:Real}
struct Bin{S, T<:Real, L<:Real}
    center::S
    width::L
    weight::T   ## Can be counts, or normalized counts or counts or normalized counts weighted by the bin length etc...
end
const Histogram{S, T, L} = Vector{Bin{S,T, L}} where {S, T, L}

function Histogram(data::AbstractVector{T}; kwargs...) where {T}
    return histogram(Data(data, 1); kwargs...)
end

function Base.extrema(H::Histogram)
    lo = minimum(getfield.(H, :center)) - first(getfield.(H, :width))/2
    hi = maximum(getfield.(H, :center)) + last(getfield.(H, :width))/2
    return lo, hi
end
function adapt_histogram(H::Histogram{S,T}; normalization=1, xscale=1) where {S<:Real, T<:Real}
    if xscale != 1
        H = map(H) do b
            Bin(b.center / xscale, b.width / xscale, b.weight)
        end
    end
    if normalization == 1
        return H 
    end
    if normalization isa Number
        H = map(H) do b
            Bin(b.center, b.width, b.weight / normalization)
        end
    elseif normalization == :probability
        total_weight = sum(b.weight for b in H)
        H = map(H) do b
            Bin(b.center, b.width, b.weight / (total_weight))
        end
    elseif normalization == :area 
        total_area = sum(b.weight * b.width for b in H)
        H = map(H) do b
            Bin(b.center, b.width, b.weight / (total_area))
        end
    else
        error("Unknown normalization type: $normalization")
    end
    return H
end

function Histogram(data::DataSet{T, L, A}; x_len = nothing, lo=nothing, hi=nothing, bins=100) where {T<:Real, L<:Real, A}
    if isnothing(lo) && isnothing(hi)
        lohi = map(data) do d
            lo, hi = ThreadsX.extrema(d.data)
            return lo, hi
        end 
        lo = minimum(first.(lohi))
        hi = maximum(last.(lohi))
    end
    if isnothing(lo)
        lo = map(data) do d
            ThreadsX.minimum(d.data)
        end |> minimum
    end
    if isnothing(hi)
        hi = map(data) do d
            ThreadsX.maximum(d.data)
        end |> maximum
    end

    isnothing(x_len) && (x_len = (hi - lo)/bins)
    return Histogram(data, lo, hi, x_len)
end

function Histogram(data::DataSet{T, L, A}, lo::Real, hi::Real, x_len::Real) where {T<:Real, L<:Real, A}
    invw = 1 / x_len
    nbins = ceil(Int,(hi - lo) * invw)
    width = x_len = (hi - lo) / nbins
    invw = 1 / width
    

    centers = Vector{Float64}(undef, nbins)
    @inbounds for i in 1:nbins
        centers[i] = lo + (i - 0.5) * width
    end

    counts = zeros(L, nbins)
    for d in data 
        Chunks = index_chunks(1:length(d.data); n = Threads.nthreads())
        local_counts = [zeros(Float64, nbins) for _ in 1:Threads.nthreads()]
        Threads.@threads for (idx,r) in enumerate(Chunks)
            @inbounds local_counts[idx] = count_appearances(@view(d.data[r]), lo, hi, invw, nbins)
        end
        local_count = reduce(+, local_counts)
        counts .+= local_count .* d.weight
    end
    total_counts = sum(length(d.data) * d.weight for d in data)
    bins = map(1:nbins) do i 
        Bin(centers[i], width, counts[i]/total_counts)
    end
    return bins
end

function count_appearances(dat::A, lo, hi, invw, bins) where {A<:AbstractVector{T}} where {T<:Real}
    local_counts = zeros(Int, bins)
    @inbounds for d in dat
        if d >= lo && d <= hi
            bin_idx = floor(Int,(d - lo) * invw) + 1
            local_counts[bin_idx] += 1
        end
    end
    return local_counts
end

function rebin(H::Histogram{S,T}; bins::Int=100, x_len=nothing) where {S<:Real, T<:Real}
    lo, hi = extrema(H)
    isnothing(x_len) && (x_len = (hi - lo)/bins)
    bins = ceil(Int,(hi - lo)/x_len)
    x_len = (hi - lo)/bins
    n = length(H)
    if bins >= n
        return H
    end
    ## Here, we assume that H describes a step function, which is constant over each bin (such that bins does not have to be a divisor of the number of bins in H)

    new_bins = Vector{Bin{S,T}}(undef, bins)
    for i in 1:bins 
        start_edge = lo + (i - 1) * x_len
        end_edge = lo + i * x_len
        total_weight = zero(T)
        total_width = zero(S)
        total_center_weight = zero(S * T)
        @inbounds for b in H
            bin_start = b.center - b.width / 2
            bin_end = b.center + b.width / 2
            overlap_start = max(start_edge, bin_start)
            overlap_end = min(end_edge, bin_end)
            if overlap_start < overlap_end
                overlap_width = overlap_end - overlap_start
                weight_fraction = overlap_width / b.width
                total_weight += b.weight * weight_fraction
                total_width += overlap_width
                total_center_weight += b.center * b.weight * weight_fraction
            end
        end
        center = total_center_weight / total_weight
        new_bins[i] = Bin(center, total_width, total_weight)
    end
    return new_bins
end
