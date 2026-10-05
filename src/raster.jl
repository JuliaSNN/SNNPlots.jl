import SNNModels: resample_spikes

"""
    raster(spiketimes::Spiketimes, t = nothing, markersize = 1)
    raster(P, t = nothing; kwargs...)

Raster plot of spike times in a new Makie figure. Returns a `Makie.FigureAxisPlot`
(destructure it as `fig, ax, plt = raster(...)`).

Two methods:

- `raster(spiketimes, t, markersize)`: `spiketimes` is a `Spiketimes`
  (`Vector{Vector{Float32}}`, one vector of spike times in ms per neuron, as returned by
  `spiketimes(pop)`). `t` is the time window `[t_start, t_end]` in ms (default
  `[0, latest spike]`); `markersize` is positional. Neuron `n` is drawn on row `n`. Spike times are
  plotted in ms, although the axis label reads "Time (s)".
- `raster(P, t; kwargs...)`: `P` is a population, a stimulus, or a `NamedTuple` of them (for
  example `model.pop`). Creates a `Figure` and an `Axis` and calls [`raster!`](@ref) with the same
  arguments; see there for the keyword arguments. Time is shown in seconds.

The populations must record spikes: call `monitor!(pop, [:fire])` before the simulation.
At most 200 000 spikes are drawn; larger sets are randomly subsampled (with a warning) by
`SNNModels.resample_spikes`.

# Example
```julia
using SpikingNeuralNetworks
SNN.@load_units
E = SNN.Poisson(N = 50, param = SNN.PoissonParameter(10Hz))
SNN.monitor!(E, [:fire])
SNN.sim!([E]; duration = 1s)
fig, ax, plt = SNN.raster(E, 0:1s)       # population method, time axis in s
fig2, ax2, plt2 = SNN.raster(SNN.spiketimes(E))  # Spiketimes method
```
"""
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
    isempty(Y) || ylims!(0, maximum(Y) + 1)
    t = typeof(t) <: AbstractRange ? t[[1, end]] : t
    return Makie.FigureAxisPlot(fig, ax, plt)
end

"""
    raster!(ax::Axis, spiketimes::Spiketimes, t = nothing; markersize = 1, order = nothing, kwargs...)
    raster!(ax, P, t = nothing; interval = nothing, populations = nothing, names = nothing,
            every = 1, markersize = 1, order = [], kwargs...)

Draw a raster plot into the existing Makie axis `ax` and return the scatter plot object.

`Spiketimes` method: `spiketimes` holds one vector of spike times (ms) per neuron; `t` is the
window `[t_start, t_end]` (ms, default `[0, latest spike]`); `order` is an optional permutation
of the neurons (neuron `order[k]` is drawn on row `k`); remaining `kwargs` are passed to
`Makie.scatter!`. Times are plotted in ms.

Population method: `P` is a population/stimulus or a `NamedTuple` of them; populations are stacked
vertically, separated by red dashed lines, and labelled on the y axis with the first 10
characters of their `name`. Keyword arguments:
- `t` or `interval`: time window in ms (a range or a two-element vector); `t` has precedence.
  The x axis is in seconds (spike times are divided by `s`).
- `populations`: if given, `P` must be a single population and `populations` a vector of index
  vectors; each index set is drawn as a separate group, labelled with `names`
  (default `"pop_i"`). Without `populations`, `names` is ignored and the populations' `name`
  fields are used.
- `every = 1`: draw one spike out of `every`.
- `markersize = 1`: marker size.
- `order = []`: neuron order passed to the per-population raster.
- other `kwargs` are ignored by the Makie backend.

Requires a `:fire` record on every population (`monitor!(pop, [:fire])`). The y-limits are set
with `ylims!` on the current axis, which is `ax` only if `ax` is the most recently created axis.

# Example
```julia
using SpikingNeuralNetworks, CairoMakie
SNN.@load_units
E = SNN.Poisson(N = 40, param = SNN.PoissonParameter(10Hz), name = "E")
I = SNN.Poisson(N = 10, param = SNN.PoissonParameter(20Hz), name = "I")
SNN.monitor!([E, I], [:fire])
SNN.sim!([E, I]; duration = 1s)
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (s)")
SNNPlots.raster!(ax, (; E, I), 0:1s)   # raster! is not exported
```
"""
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
    isempty(Y) || ylims!(ax, 0, maximum(Y) + 1)
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
    interval = nothing,
    populations = nothing,
    names = nothing,
    every = 1,
    markersize = 1,
    order::Vector = [],
    kwargs...,
)
    t = isnothing(t) ? interval : t
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
        isempty(Y) || ylims!(0, maximum(Y) + 1)
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
    t = isnothing(t) ? nothing : t[[1, end]]
    X, Y = Float32[], Float32[]
    order = isempty(order) ? eachindex(spiketimes) : order
    for n in order
        for st in spiketimes[n]
            if isnothing(t) || (st > t[1] && st < t[2])
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
