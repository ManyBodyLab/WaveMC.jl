struct Data{T, L<:Real, A<:AbstractVector{T}}
    data::A
    weight::L
end
const DataSet{T, L, A} = Vector{Data{T,L,A}} where {A<:AbstractVector{T}} where {T, L<:Real}

linear_edges(lo, hi, n::Int) = collect(range(lo, hi, length=n+1))


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
    density::Bool = false,
    xscale::Real = 1,
) where {T<:Real, L<:Real, A}
    density && (normalization = :pdf)
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

    data_all = reduce(vcat, d.data for d in dataset)
    wts_all  = reduce(vcat, fill(Float64(d.weight), length(d.data)) for d in dataset)
    h = fit(StatsBase.Histogram, data_all, StatsBase.weights(wts_all), edge_vec)
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
    new_weights = zeros(Float64, length(edge_vec) - 1)
    tmp = StatsBase.Histogram(edge_vec)
    for (c, w) in zip(centers, float.(h.weights))
        idx = StatsBase.binindex(tmp, c)
        checkbounds(Bool, new_weights, idx) && (new_weights[idx] += w)
    end
    return StatsBase.Histogram(edge_vec, new_weights)
end

function histogram_interpolation(
    h::StatsBase.Histogram;
    order::Int = 1,
    extrapolation_bc = Interpolations.Flat()
)
    centers = bin_centers(h)
    density = bin_weights(h) ./ bin_widths(h)
    if order == 0
        itp = Interpolations.interpolate((centers,), density, Interpolations.Gridded(Interpolations.Constant()))
        return Interpolations.extrapolate(itp, extrapolation_bc)
    elseif order == 1
        return Interpolations.linear_interpolation(centers, density; extrapolation_bc)
    else
        error("order must be 0 or 1")
    end
end
