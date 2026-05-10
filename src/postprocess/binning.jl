struct Data{T, L<:Real, A<:AbstractVector{T}}
    data::A
    weight::L
end
const DataSet{T, L, A} = Vector{Data{T,L,A}} where {A<:AbstractVector{T}} where {T, L<:Real}

linear_edges(lo, hi, n::Int) = collect(range(lo, hi, length=n+1))

function Histogram(
    dataset::DataSet{T,L,A};
    bins::Int = ceil(Int, log2(length(dataset))) + 1,
    lo = nothing,
    hi = nothing,
    edges = nothing,
    normalization = :none,
    density::Bool = false,
    total = 1,
    xscale::Real = 1,
    area = v->diff(v),
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
    areas = area(edge_vec)

    nbins = length(edge_vec) - 1
    lo_edge = first(edge_vec)
    hi_edge = last(edge_vec)
    weights = extract_data_weights(dataset, edge_vec, nbins, lo_edge, hi_edge)
    l = sum(weights)
    weights = (weights ./ l).* total ./areas
    h = StatsBase.Histogram(edge_vec, weights)
    return normalize(h, mode=normalization)
end
function extract_data_weights(dataset::DataSet{T,L,A}, edge_vec, nbins, lo_edge, hi_edge) where {T<:Real, L<:Real, A}
    weights = zeros(L, nbins)
    for data in dataset
        weights .+= extract_weights(data, edge_vec, nbins, lo_edge, hi_edge)
    end
    return weights
end
function extract_weights(data::Data{T,L,A}, edge_vec, nbins, lo_edge, hi_edge) where {T<:Real, L<:Real, A}
    step_ = edge_vec isa AbstractRange ? step(edge_vec) : 0.0
    chunks = index_chunks(eachindex(data.data); n = Threads.nthreads())
    thread_weights = [zeros(Int, nbins) for _ in 1:Threads.nthreads()]
    Threads.@threads for (idx, r) in enumerate(chunks)
        thread_weights[idx] = _extract_weights(data.data, edge_vec, nbins, lo_edge, hi_edge, step_, r)
    end
    return sum(thread_weights) * data.weight
end
function _extract_weights(data, edge_vec, nbins, lo_edge, hi_edge, step_, indices)
    local_weights = zeros(Int, nbins)
    @inbounds @simd for i in indices
        x = data[i]
        if lo_edge <= x <= hi_edge
            idx = _find_bin(edge_vec, x, nbins, lo_edge, step_)
            local_weights[idx] += 1
        end
    end
    return local_weights
end
@inline function _find_bin(edge_vec, x, nbins, lo_edge, step_)
    return clamp(searchsortedlast(edge_vec, x), 1, nbins)
end
@inline function _find_bin(edge_vec::AbstractRange, x, nbins, lo_edge, step_)
    return clamp(floor(Int, (x - lo_edge) / step_) + 1, 1, nbins)
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
    lo_edge = first(edge_vec)
    hi_edge = last(edge_vec)
    for (c, w) in zip(centers, float.(h.weights))
        c < lo_edge || c > hi_edge && continue
        idx = clamp(searchsortedlast(edge_vec, c), 1, length(new_weights))
        new_weights[idx] += w
    end
    return StatsBase.Histogram(edge_vec, new_weights)
end

function histogram_interpolation(
    h::StatsBase.Histogram;
    order::Int = 1,
    extrapolation_bc = Interpolations.Flat(),
    density = false,
)
    centers = bin_centers(h)
    dens = bin_weights(h)
    if density 
        widths = bin_widths(h)
        dens ./= widths
    end
    knots = centers

    if order == 0
        edge_knots = collect(first(h.edges))
        dens_ext = [first(dens); dens]
        return extrapolate(interpolate((edge_knots,), dens_ext, Gridded(Constant{Previous}())), extrapolation_bc)
    elseif order == 1
        return extrapolate(interpolate((collect(knots),), dens, Gridded(Linear())), extrapolation_bc)
        ## TODO: Use order = 1 extrapolation to use regular grid for higher order interpolation
    # elseif order == 2
    #     return extrapolate(Interpolations.scale(interpolate(dens, BSpline(Quadratic(Flat(OnGrid())))), (knots,)), extrapolation_bc)
    # elseif order == 3
    #     return cubic_spline_interpolation(centers, dens; extrapolation_bc=extrapolation_bc)
    #     return extrapolate(Interpolations.scale(interpolate(dens, BSpline(Cubic(Flat(OnGrid())))), (knots,)), extrapolation_bc)
    else
        error("order must be 0, 1")
    end
end
