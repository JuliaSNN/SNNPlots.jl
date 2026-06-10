import SNNModels: resample_spikes

function raster(spiketimes::Spiketimes, t = nothing, markersize=1)
    t = isnothing(t) ? [0, maximum(vcat(spiketimes...))] : t
    X, Y = _raster(spiketimes, t)
    X, Y = resample_spikes(X, Y)
    fig, ax, plt = scatter(
        X,
        Y,
        markersize = markersize,
        color = :black,
        axis = (;
            xlabel = "Time (s)",
            ylabel = "Neuron",
        ),
    )
    xlims!(ax, extrema(t))
    ylims!(0, maximum(Y) + 1)
    t = typeof(t) <: AbstractRange ? t[[1, end]] : t
    return Makie.FigureAxisPlot(fig, ax, plt)
end

function raster!(ax::Axis, spiketimes::Spiketimes, t = nothing; markersize=1, order=nothing, kwargs...)
    t = isnothing(t) ? [0, maximum(vcat(spiketimes...))] : t
    order = isnothing(order) ? eachindex(spiketimes) : order

    X, Y = _raster(spiketimes[order], t)
    X, Y = resample_spikes(X, Y)
    plt = scatter!(
        ax,
        X,
        Y;
        markersize = markersize,
        color = :black,
        kwargs...,
    )
    t = typeof(t) <: AbstractRange ? t[[1, end]] : t
    !isnothing(t) && xlims!(ax, extrema(t))
    ylims!(ax, 0, maximum(Y) + 1)
    return plt
end

function raster(P, t = nothing; kwargs...)
    if _backend == :Plots
        ax = plot(
            m = (1, :black),
            leg = :none,
            xaxis = ("Time (s)", (0, Inf)),
            yaxis = ("Neuron",),
            label = "",
        )
        return raster!(ax, P, t; kwargs...)
    elseif _backend == :Makie
        f = Figure()
        ax = Axis(f[1, 1], 
            xlabel = "Time (s)",
            ylabel = "Neuron",
        )
        plt = raster!(ax, P, t; kwargs...)
        return Makie.FigureAxisPlot(f, ax, plt)
    end
end


function raster!(
    ax,
    P,
    t = nothing;
    populations = nothing,
    names = nothing,
    every = 1,
    markersize = 1,
    order::Vector = [],
    kwargs...,
)
    if isnothing(populations)
        y0 = Int32[0]
        X = Float32[]
        Y = Float32[]
        names = Vector{String}()
        P = typeof(P) <: AbstractPopulation ? [P] : [getfield(P, k) for k in keys(P)]
        for p in P
            x, y, _y0 = _raster(p, t; order)
            push!(names, p.name)
            append!(X, x)
            append!(Y, y .+ sum(y0))
            isempty(_y0) ? push!(y0, p.N) : (y0 = vcat(y0, _y0))
        end
    else
        @assert typeof(P) <: AbstractPopulation
        X, Y, y0 = _raster_populations(P, t; populations = populations)
        names = isnothing(names) ? ["pop_$i" for i = 1:length(P)] : names
    end

    X, Y = resample_spikes(X, Y)
    X = X ./ s

    plt = scatter!(
        ax,
        X[1:every:end],
        Y[1:every:end];
        color = :black,
        markersize
    )
    t = typeof(t) <: AbstractRange ? t[[1, end]] : t
    if _backend == :Plots
        !isnothing(t) && plot!(xlims = t ./ s)
        plot!(yticks = (cumsum(y0)[1:(end-1)] .+ (y0 ./ 2)[2:end], names), yrotation = 45)
        y0 = y0[2:(end-1)]
        !isempty(y0) && hline!(ax, cumsum(y0), linecolor = :red, label = "")
        plot!(ax; kwargs...)
        return ax
    elseif _backend == :Makie
        !isnothing(t) && Makie.xlims!(ax, t ./ s)
        ax.yticks = (cumsum(y0)[1:(end-1)] .+ (y0 ./ 2)[2:end], [n[1:minimum([10, length(n)])] for n in string.(names)])
        ax
        y0 = y0[1:(end)]
        !isempty(y0) && Makie.hlines!(ax, cumsum(y0), color = :red, linewidth = 1, label = "", linestyle = :dash)
        ylims!(0, maximum(Y) + 1)
        return plt
    end
end

function _raster_populations(
    p,
    t = nothing;
    populations::Vector{T},
) where {T<:AbstractVector}
    all_spiketimes = spiketimes(p)
    x, y = Float32[], Float32[]
    y0 = Int32[0]
    for pop in populations
        spiketimes_pop = all_spiketimes[pop] ## population spiketimes
        for n in eachindex(spiketimes_pop) ## neuron spiketimes
            for st in spiketimes_pop[n] ## spiketime
                if isnothing(st) || (st > t[1] && st < t[end])
                    push!(x, st)
                    push!(y, n + cumsum(y0)[end])
                end
            end
        end
        push!(y0, length(spiketimes_pop))
    end
    return x, y, y0
end

function _raster(spiketimes::Spiketimes, t = nothing; order = [])
    t, X, Y = t[[1, end]], Float32[], Float32[]
    order = isempty(order) ? eachindex(spiketimes) : order
    for n in order
        for st in spiketimes[n]
            if isnothing(st) || (st > t[1] && st < t[2])
                push!(X, st)
                push!(Y, n)
            end
        end
    end
    return X, Y
end





function _raster(
    p::T,
    interval = nothing;
    order = [],
) where {T<:Union{AbstractPopulation,AbstractStimulus}}
    !haskey(p.records, :fire) && @error "No fire record found in population $(p.name)"
    interval = typeof(interval) <: AbstractRange ? interval[[1, end]] : interval
    st = SNNModels.spiketimes(p; interval)
    x, y = _raster(st, interval; order)
    return x, y, Int32[]
end
