function Histogram(filedir, observable;
        bins = 100,
        normalization=1,
        xscale=1,
        x_len=nothing,
        kwargs...
    )
    tasks = parse.(Int,replace.(filter(x->startswith(x, "task"), readdir(filedir)), "task" => ""))

    !(normalization isa Vector) && (normalization = fill(normalization, length(tasks)))
    !(xscale isa Vector) && (xscale = fill(xscale, length(tasks)))
    !(x_len isa Vector) && (x_len = fill(x_len, length(tasks)))
    !(bins isa Vector) && (bins = fill(bins, length(tasks)))

    histograms = map(eachindex(tasks)) do i
        Histogram(filedir, observable, tasks[i]; normalization=normalization[i], xscale=xscale[i], x_len=x_len[i], bins=bins[i], kwargs...)
    end
    return histograms
end
function Histogram(filedir, observable::AbstractVector{<:AbstractString}, task::Int;
    x_len=nothing,
    xscale=1,
    bins = 100,
    lo=nothing,
    hi=nothing,
    kwargs...
)
    taskname = "task$(lpad(string(task), 4, '0'))"
    file = only(filter(x->x == taskname, readdir(filedir)))

    runs = filter!(x->startswith(x, "run") && occursin("meas",x),readdir(joinpath(filedir,file)))

    x = [read_runs(filedir, taskname, runs, obs) for obs in observable]
    l = length(x[1])
    bin_lengths = only(unique([d[2] for d in x])) ## Bin lengths of observables should be the same!
    L = length(observable)
    
    dataset = map(1:l) do i 
        datas = x[i]
        weight = bin_lengths[i]
    end


    datas = Iterators.product((d[1] for d in x)...)
    


    dataset = [Data(datas[i], bin_lengths[i]) for i in eachindex(datas)]
    x_len = isnothing(x_len) ? nothing : x_len*xscale
    return adapt_histogram(Histogram(dataset; x_len=x_len, lo=lo, hi=hi, bins=bins); xscale=xscale, kwargs...)
end

function Histogram(filedir, observable::AbstractString, task::Int;
    x_len=nothing,
    xscale=1,
    bins = 100,
    lo=nothing,
    hi=nothing,
    kwargs...
)
    taskname = "task$(lpad(string(task), 4, '0'))"
    file = only(filter(x->x == taskname, readdir(filedir)))

    runs = filter!(x->startswith(x, "run") && occursin("meas",x),readdir(joinpath(filedir,file)))

    datas, bin_lengths = read_runs(filedir, taskname, runs, observable)
    dataset = [Data(datas[i], bin_lengths[i]) for i in eachindex(datas)]
    x_len = isnothing(x_len) ? nothing : x_len*xscale
    return adapt_histogram(Histogram(dataset; x_len=x_len, lo=lo, hi=hi, bins=bins); xscale=xscale, kwargs...)
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

function plot_histogram(H::Histogram;
    xlabel=L"r",
    ylabel=L"g(r)",
    rescaling_function = x->1, ## Rescales the y-value by rescaling_function(center)
    title="",
    color=:blue,
)
    centers = [b.center for b in H]
    data = [b.weight/b.width for b in H]
    data = rescaling_function.(centers) .* data
    lo, hi = extrema(H)

    ylim = maximum(data)

    p = sortperm(centers)
    data, centers = data[p], centers[p]

    fig = with_theme(theme_latexfonts()) do 
        fig = Figure(size = (800, 600),fontsize=20)
        ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, title = title)
        CairoMakie.scatter!(ax, centers, data; color=color,strokecolor=color)
        xlims!(ax, lo, hi)
        ylims!(ax, 0, ylim * 1.1)
        return fig
    end
    return fig
end
