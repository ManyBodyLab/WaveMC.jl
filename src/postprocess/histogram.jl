function Histogram(filedir, observable;
        bins = 100,
        normalization = :none,
        xscale = 1,
        edges = nothing,
        kwargs...
    )
    tasks = parse.(Int, replace.(filter(x->startswith(x, "task"), readdir(filedir)), "task" => ""))

    !(normalization isa Vector) && (normalization = fill(normalization, length(tasks)))
    !(xscale isa Vector) && (xscale = fill(xscale, length(tasks)))
    !(bins isa Vector) && (bins = fill(bins, length(tasks)))
    !(edges isa Vector) && (edges = fill(edges, length(tasks)))

    histograms = map(eachindex(tasks)) do i
        Histogram(filedir, observable, tasks[i]; normalization=normalization[i], xscale=xscale[i], bins=bins[i], edges=edges[i], kwargs...)
    end
    return histograms
end

function Histogram(filedir, observable::AbstractVector{<:AbstractString}, task::Int;kwargs...)
    taskname = "task$(lpad(string(task), 4, '0'))"
    file = only(filter(x->x == taskname, readdir(filedir)))
    runs = filter!(x->startswith(x, "run") && occursin("meas",x), readdir(joinpath(filedir,file)))

    x = [read_runs(filedir, taskname, runs, obs) for obs in observable]
    bin_lengths = only(unique([d[2] for d in x]))

    datas = Iterators.product((d[1] for d in x)...)
    dataset = [Data(collect(datas[i]), bin_lengths[i]) for i in eachindex(datas)]
    return Histogram(dataset; kwargs...)
end

function Histogram(filedir, observable::AbstractString, task::Int; kwargs...)
    taskname = "task$(lpad(string(task), 4, '0'))"
    file = only(filter(x->x == taskname, readdir(filedir)))
    runs = filter!(x->startswith(x, "run") && occursin("meas",x), readdir(joinpath(filedir,file)))

    datas, bin_lengths = read_runs(filedir, taskname, runs, observable)
    dataset = [Data(datas[i], bin_lengths[i]) for i in eachindex(datas)]
    return Histogram(dataset; kwargs...)
end

function read_runs(filedir, taskname, runs, observable)
    data, bin_length = h5open(joinpath(filedir,taskname,runs[1]), "r") do fid
        dat = read(fid["observables"][observable]["samples"])
        bin_len = read(fid["observables"][observable]["bin_length"])
        return dat, bin_len
    end
    datas = Vector{Vector{eltype(data)}}(undef, length(runs))
    datas[1] = data
    bin_lengths = [bin_length]
    for r in eachindex(runs)[2:end]
        h5open(joinpath(filedir,taskname,runs[r]), "r") do fid
            datas[r] = read(fid["observables"][observable]["samples"])
            push!(bin_lengths, read(fid["observables"][observable]["bin_length"]))
        end
    end
    return datas, bin_lengths
end

function plot_histogram(h::StatsBase.Histogram;
    xlabel = L"r",
    ylabel = L"g(r)",
    rescaling_function = x -> 1,
    title = "",
    color = :blue,
    density = false,
    area_function = bin_widths,
)
    @show rescaling_function
    centers = bin_centers(h)
    ydata = density ? bin_weights(h) ./ area_function(h) : float.(bin_weights(h))
    ydata = rescaling_function.(centers) .* ydata

    lo = first(first(h.edges))
    hi = last(first(h.edges))

    p = sortperm(centers)

    fig = with_theme(theme_latexfonts()) do
        fig = Figure(size=(800, 600), fontsize=20)
        ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title)
        CairoMakie.scatter!(ax, centers[p], ydata[p]; color=color, strokecolor=color)
        xlims!(ax, lo, hi)
        ylims!(ax, 0, maximum(ydata) * 1.1)
        return fig
    end
    return fig
end
