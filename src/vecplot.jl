# Vector plot
using Statistics: quantile
"""
    vecplot(p, sym::Symbol; kwargs...)
    vecplot(p, syms::Vector{Symbol}; interval = nothing, ylabel = "", title = "", kwargs...)
    vecplot(p, sym, interval::AbstractRange; kwargs...)

Plot the recorded time course of the variable(s) `sym` of the population (or connection,
stimulus) `p` in a new Makie figure; returns a `Makie.FigureAxisPlot`.

The variable must have been recorded with `monitor!(p, [sym])` before the simulation. One line is
drawn per neuron (or one line for the population mean with `pop_average = true`); every keyword
argument not listed here is forwarded to [`vecplot!`](@ref) (`neurons`, `pop_average`, `r`,
`sym_id`, `factor`, `add_spikes`, `variables`, `lw`, `color`, `ribbon`, ...).

- `interval`: time range in ms (e.g. `0:1:1000ms` or `0:1s`); default: the whole record.
- `ylabel`, `title`: axis labels.

The x axis data are in ms and the tick labels are shown in seconds.

The methods `vecplot(P::Array, sym)` and `vecplot(P, syms::Array)` (one panel per population or
per symbol) call the Plots.jl `plot(...; layout)` API and do not work with the Makie backend used
by SNNPlots 0.2.10.

# Example
```julia
using SpikingNeuralNetworks
SNN.@load_units
E = SNN.IF(N = 5, param = SNN.IFParameter(gl = 10nS, C = 200pF))
SNN.monitor!(E, [:v])
stim = SNN.CurrentStimulus(E; param = SNN.CurrentNoise(E; I_base = 300pA))
model = SNN.compose(; E, stim)
SNN.sim!(; model, duration = 200ms)
fig, ax, plt = SNN.vecplot(E, :v, interval = 0:1:200ms, neurons = 1:2, ylabel = "V (mV)")
```
"""
function vecplot(p, sym::Symbol; kwargs...)
    vecplot( p, [sym]; kwargs...)
end

function vecplot(p, sym::Vector{Symbol}; interval=nothing, ylabel="", title="", kwargs...)
    if _backend == :Plots
        ax = plot(; xlabel = "Time (s)", ylabel = ylabel, title = title)
        for s in sym
            vecplot!(ax, p, s; interval, label = string(s), kwargs...)
        end
        plot!(ax, ylims = :auto)
        return ax
    elseif _backend == :Makie
        f = Figure()
        ax = Axis(f[1, 1],
            xlabel = "Time (s)";
            ylabel ,
            title 
        )
        plt = nothing
        for s in sym
            plt = vecplot!(ax, p, s; interval, labels = [string(s)], kwargs...)
        end
        return Makie.FigureAxisPlot(f, ax, plt)
    end
end

function vecplot(P::Array, sym; kwargs...)
    plts = [vecplot(p, sym; kwargs...) for p in P]
    N = length(plts)
    plot(plts..., size = (600, 400N), layout = (N, 1))
end

vecplot(p, sym, interval::T; kwargs...) where {T<:AbstractRange} =
    vecplot(p, sym; interval = interval, kwargs...)

vecplot!(my_plot, p, sym::Symbol, interval::T; kwargs...) where {T<:AbstractRange} =
    vecplot!(my_plot, p, sym; interval = interval, kwargs...)



function _match_r(r, r_v)
    r = isnothing(r) ? range(r_v[1], r_v[end]) : r
    r[end] > r_v[end] &&
        throw(ArgumentError("The end time is greater than the record time"))
    r[1] < r_v[1] && throw(ArgumentError("The start time is less than the record time"))
    return r
end

"""
    vecplot!(ax, p, sym; neurons = nothing, pop_average = false, interval = nothing, r = nothing,
             sym_id = nothing, factor = 1.0f0, add_spikes = false, variables = nothing, lw = 2,
             color = nothing, ribbon = false, label = nothing, kwargs...)
    vecplot!(ax, p, sym::Symbol, interval::AbstractRange; kwargs...)

Draw the recorded variable `sym` of `p` into the Makie axis `ax`; returns the last `lines!` plot.

The record is read with `SNNModels.record(p, sym; variables, range = true)`, which returns an
interpolated record and its time range.

Keyword arguments:
- `interval` (or `r`; `interval` has precedence): time range in ms at which the interpolated record
  is sampled. Default: from the first to the last recorded time in steps of 1 ms. An error is
  thrown if the range exceeds the recorded time.
- `neurons`: neuron indices (an `Int` or a collection); default all.
- `pop_average = false`: plot the mean over the selected neurons.
- `ribbon = false`: with `pop_average = true`, also draw a band between the 20th and 80th
  percentiles across neurons.
- `sym_id`: required when the record is three-dimensional (e.g. one value per dendritic
  compartment): index of the second dimension to plot.
- `factor = 1.0f0`: scalar multiplying the trace.
- `add_spikes = false`: set the trace to 20 mV at the first sample after every spike (needs a
  `:fire` record).
- `variables`: forwarded to `record` (selects variables of a nested record).
- `lw = 2`: line width; `color`: line colour (default: cycle through the theme palette by neuron
  index); `label`: legend label (converted with `string`).
- other `kwargs` are ignored by the Makie backend.

x data are in ms; the six x ticks are relabelled in seconds.

# Example
```julia
using SpikingNeuralNetworks, CairoMakie
SNN.@load_units
E = SNN.IF(N = 5, param = SNN.IFParameter(gl = 10nS, C = 200pF))
SNN.monitor!(E, [:v])
stim = SNN.CurrentStimulus(E; param = SNN.CurrentNoise(E; I_base = 300pA))
model = SNN.compose(; E, stim)
SNN.sim!(; model, duration = 200ms)
fig = Figure()
ax = Axis(fig[1, 1], ylabel = "V (mV)")
SNN.vecplot!(ax, E, :v; pop_average = true, ribbon = true)
```
"""
function vecplot!(
    ax,
    p,
    sym;
    neurons = nothing,
    pop_average = false,
    interval = nothing,
    r = nothing,
    sym_id = nothing,
    factor = 1.0f0,
    add_spikes = false,
    variables = nothing,
    lw = 2,
    color = nothing,
    ribbon = false,
    label = nothing,
    kwargs...,

)
    # get the record and its sampling rate
    y, r_v = SNNModels.record(p, sym; variables, range = true)
    r = isnothing(interval) ? r : interval
    r = _match_r(r, r_v)

    neurons = isnothing(neurons) ? axes(y, 1) : neurons
    neurons = isa(neurons, Int) ? [neurons] : neurons
    
    if isa(factor, Symbol)
        factor, _ = SNN.interpolated_record(p, factor)
        factor = factor(neurons, r)
    elseif isa(factor, Matrix)
        factor = factor(neurons, :)
        @assert size(factor, 1) == length(neurons) "The factor matrix must have the same number of rows as the number of neurons"
        @assert size(factor, 2) == size(y, 2) "The factor matrix must have the same number of columns as the number of time points in the record"
    end


    # check if the record is a vector or a matrix
    if ndims(y) == 3
        isnothing(sym_id) && (throw(
            ArgumentError(
                "The record is a matrix, please specify the index ($sym_id) of the matrix to plot with `sym_id`",
            ),
        ))
        y = y(neurons, sym_id, r)
    else
        y = y(neurons, r)
    end
    

    
    band_up = pop_average && ribbon ? [quantile(y[:, t], 0.8) for t in axes(y, 2)] : nothing
    band_down = pop_average && ribbon ? [quantile(y[:, t], 0.2) for t in axes(y, 2)] : nothing
    y = pop_average ? SNNModels.Statistics.mean(y, dims = 1) : y




    if add_spikes
        @assert haskey(p.records, :fire) "No fire record found in population $(p.name)"
        spiketimes = SNNModels.spiketimes(p)[neurons]
        @assert length(spiketimes) == size(y, 1) "The number of spiketimes $(length(spiketimes)) does not match the number of neurons $(size(y, 1)) in population $(p.name)"
        for n in eachindex(spiketimes)
            for sp in spiketimes[n]
                tt = findfirst(r .>= sp)
                isnothing(tt) && continue
                y[n, tt] = 20mV
            end
        end
    end

    @debug "Vector plot in: $(r[1])ms to $(round(Int, r[end]))ms"
    if _backend == :Plots
        return plot!(
            ax,
            r ./ 1000,
            y' .* factor,
            ribbon = ribbon,
            leg = :none,
            xaxis = ("t", extrema(r ./ 1000)),
            yaxis = (string(sym), extrema(y));
            lw = 3,
            kwargs...,
        )
        return plot!(; kwargs...)
    elseif _backend == :Makie
        my_color = isnothing(color) ? x->Cycled(x) : x->color
        isnothing(band_down) || band!(
            ax,
            r,
            (band_down) .* factor',
            (band_up) .* factor',
            color = my_color(1),
            alpha= 0.7,
        )
        plt = nothing
        for n in axes(y, 1)
            plt =  lines!(
                ax,
                r,
                y[n, :] .* factor',
                color = my_color(n),
                linewidth = lw,
                label = string(label),
            )
        end

        ax.xticks = (range(extrema(r)...,6), string.(round.(range((extrema(r)./1000)..., 6), digits=2)))
        return plt
    end
end

function vecplot(P, syms::Array; kwargs...)
    plts = [vecplot(P, sym; kwargs...) for sym in syms]
    N = length(plts)
    plot(plts..., size = (600, 400N), layout = (N, 1))
end


function raster_firing(model; path = nothing, τ = 20ms, every = 1)
    fr_average, r, labels =
        firing_rate(model.pop, interval = 0s:get_time(model), τ = τ, pop_average = true)
    pr = raster(
        model.pop,
        (get_time(model)-5s):get_time(model),
        every = every,
        size = (1200, 1000),
    )
    pf = plot(
        r,
        fr_average,
        labels = hcat(labels...),
        xlabel = "Time (s)",
        ylabel = "Firing rate (Hz)",
    )
    presults = plot(pr, pf, layout = (2, 1))
    if !isnothing(path)
        savefig(presults, path)
    end
end


export vecplot, vecplot!

