struct Data{T, L<:Real, A<:AbstractVector{T}}
    data::A
    weight::L
end
const DataSet{T, L, A} = Vector{Data{T,L,A}} where {A<:AbstractVector{T}} where {T, L<:Real}

linear_edges(lo, hi, n::Int) = collect(range(lo, hi, length=n+1))

function _fit_histogram(dataset::DataSet{T,L,A}, edges::AbstractVector; closed::Symbol=:left) where {T<:Real, L<:Real, A}
    n_total = sum(length(d.data) for d in dataset)
    data_all = Vector{T}(undef, n_total)
    wts_all  = Vector{Float64}(undef, n_total)
    idx = 1
    for d in dataset
        n = length(d.data)
        data_all[idx:idx+n-1] .= d.data
        wts_all[idx:idx+n-1]  .= d.weight
        idx += n
    end
    return fit(StatsBase.Histogram, data_all, StatsBase.weights(wts_all), edges; closed)
end

function _apply_normalization(h::StatsBase.Histogram, normalization)
    normalization in (1, :none) && return h
    normalization isa Symbol && return StatsBase.normalize(h; mode=normalization)
    h2 = float(h)
    h2.weights ./= normalization
    return h2
end

function Histogram(
    dataset::DataSet{T,L,A};
    bins::Int = 100,
    lo = nothing,
    hi = nothing,
    edges = nothing,
    normalization = :none,
    xscale::Real = 1,
) where {T<:Real, L<:Real, A}
    if isnothing(lo) && isnothing(hi)
        lohi = map(d -> ThreadsX.extrema(d.data), dataset)
        lo = minimum(first.(lohi))
        hi = maximum(last.(lohi))
    elseif isnothing(lo)
        lo = minimum(ThreadsX.minimum(d.data) for d in dataset)
    elseif isnothing(hi)
        hi = maximum(ThreadsX.maximum(d.data) for d in dataset)
    end

    lo /= xscale
    hi /= xscale

    edge_vec = if isnothing(edges)
        linear_edges(lo, hi, bins)
    elseif edges isa Function
        edges(lo, hi, bins)
    else
        collect(edges)
    end

    h = _fit_histogram(dataset, edge_vec)
    return _apply_normalization(h, normalization)
end

bin_centers(h::StatsBase.Histogram) = StatsBase.midpoints(first(h.edges))
bin_widths(h::StatsBase.Histogram)  = diff(collect(first(h.edges)))
bin_weights(h::StatsBase.Histogram) = h.weights

function rebin(h::StatsBase.Histogram; bins::Int=100, edges=nothing)
    lo = first(first(h.edges))
    hi = last(first(h.edges))
    edge_vec = if isnothing(edges)
        linear_edges(lo, hi, bins)
    elseif edges isa Function
        edges(lo, hi, bins)
    else
        collect(edges)
    end
    centers = bin_centers(h)
    wts = StatsBase.weights(float.(h.weights))
    return fit(StatsBase.Histogram, centers, wts, edge_vec)
end

function histogram_interpolation(
    h::StatsBase.Histogram;
    extrapolation_bc = Interpolations.Flat()
)
    centers = bin_centers(h)
    density = bin_weights(h) ./ bin_widths(h)
    return Interpolations.linear_interpolation(centers, density; extrapolation_bc)
end
