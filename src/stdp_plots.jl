"""
    stdp_test(stdp_param; ΔT)

Measure the weight change produced by the long-term plasticity rule `stdp_param` (any
`LTPParameter`, e.g. `STDPGerstner()`) for a single pre/post spike pair with time difference
`ΔT = t_post - t_pre` (ms).

Builds a presynaptic and a postsynaptic `Identity` neuron, each driven by its own
`SpikeTimeStimulusIdentity`: the presynaptic neuron spikes at 200 ms and the postsynaptic one
at `200ms + ΔT`. A `SpikingSynapse` pre -> post with initial weight 1 and
`LTPParam = stdp_param` is trained with `train!` for 400 ms at `dt = 0.1ms` (plasticity only
runs under `train!`). The synapse transmits into a dummy buffer, so it does not make the
postsynaptic neuron fire: exactly one pre/post pair is measured. Returns the weight change
`W - 1`.

!!! note "Changed after SNNPlots 0.2.10"
    Up to SNNPlots 0.2.10 both neurons belonged to one `Identity` population and the measured
    synapse also drove the postsynaptic neuron, which emitted an extra spike one step after
    the presynaptic spike: every measured change contained an extra causal pair at
    `Δt ≈ 0.1 ms` (for `STDPGerstner()`: `+1.6e-4` at `ΔT = 10ms` and `+3.9e-5`
    instead of a depression at `ΔT = -10ms`).

# Example
```julia
using SpikingNeuralNetworks
SNN.@load_units
dw = SNNPlots.stdp_test(SNN.STDPGerstner(); ΔT = 10ms)
```
"""
function stdp_test(stdp_param; ΔT)
    pre = Identity(N = 1, name = "pre")
    post = Identity(N = 1, name = "post")
    stim_pre = SpikeTimeStimulusIdentity(pre, :g, param = SpikeTimeParameter([200ms], [1]))
    stim_post = SpikeTimeStimulusIdentity(post, :g, param = SpikeTimeParameter([200ms + ΔT], [1]))
    w = ones(Float32, 1, 1)
    syn = SpikingSynapse(pre, post, :g, conn = w, LTPParam = stdp_param)
    syn.g = zeros(Float32, post.N) # transmit into a dummy buffer: post fires only from its stimulus
    model = compose(; pre, post, stim_pre, stim_post, syn, silent = true)
    train!(model = model, duration = 400ms, dt = 0.1ms)
    return model.syn[1].W[1] - 1
end

export stdp_test


"""
    stdp_kernel!(ax, stdp_param; ΔTs = vcat(range(-200, -0.1, 20), range(0.1, 200, 20)), fill = false, kwargs...)

Plot the STDP kernel of the plasticity rule `stdp_param` into the Makie axis `ax`.

For each `ΔT` in `ΔTs` (ms, `t_post - t_pre`) the weight change is measured with
[`stdp_test`](@ref) (a full two-neuron `train!` run per point, run in parallel with
`Threads.@threads`). The negative and positive branches are drawn as two lines with a band down to
zero. `fill` and `kwargs` are accepted but unused. Returns the last `band!` plot.

# Arguments
- `ax`: a Makie `Axis`.
- `stdp_param`: a long-term plasticity parameter (`LTPParameter`), e.g. `STDPGerstner()`.
- `ΔTs`: spike-time differences in ms (default 20 points in [-200, -0.1] and 20 in [0.1, 200]).

# Example
```julia
using SpikingNeuralNetworks, CairoMakie
SNN.@load_units
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Δt (ms)", ylabel = "ΔW")
SNN.stdp_kernel!(ax, SNN.STDPGerstner(); ΔTs = [-20.0, -5.0, 5.0, 20.0])
```
"""
function stdp_kernel!(
    ax,
    stdp_param;
    ΔTs = vcat(range(-200, -0.1, 20), range(0.1, 200, 20)),
    fill = false,
    kwargs...,
)
    ΔWs = zeros(Float32, length(ΔTs))
    Threads.@threads for i in eachindex(ΔTs)
        ΔT = ΔTs[i]
        ΔWs[i] = stdp_test(stdp_param; ΔT)
    end

    n_plus = findall(ΔTs .>= 0)
    n_minus = findall(ΔTs .< 0)
    lines!(
        ax,
        ΔTs[n_minus],
        ΔWs[n_minus],
    )
    lines!(
        ax,
        ΔTs[n_plus],
        ΔWs[n_plus],
    )
    band!(ax, ΔTs[n_minus], zeros(size(n_minus)),ΔWs[n_minus])
    band!(ax, ΔTs[n_plus], zeros(size(n_plus)), ΔWs[n_plus])
    # plot!(ylims = extrema(ΔWs) .* 1.4, xlims = extrema(ΔTs), framestyle = :zerolines)
    # plot!(; kwargs...)
    # plot!(ylims = :auto, xlims = extrema(ΔTs))
end

"""
    stdp_kernel(stdp_param; kwargs...)

Create a 500x300 Makie `Figure` with an axis titled "STDP" (x: `T_post - T_pre`, y: `ΔW`), draw
the kernel of `stdp_param` with [`stdp_kernel!`](@ref) (same keyword arguments) and return the
`Figure`.

# Example
```julia
using SpikingNeuralNetworks
SNN.@load_units
fig = SNN.stdp_kernel(SNN.STDPGerstner(); ΔTs = [-20.0, -5.0, 5.0, 20.0])
```
"""
function stdp_kernel(args...; kwargs...) 
    fig = Figure(size = (500, 300))
    ax = Axis(fig[1, 1], 
        xlabel = "T_post - T_pre ",
        ylabel = "ΔW",
        title = "STDP",
    )
    stdp_kernel!(ax, args...; kwargs...)
    return fig
end

# function stdp_weight_correlated(stdp_param, rate1, rate2, τ_cov = 10ms)
#     T = 120_000ms
#     N_spike1 = T * rate1 |> Int
#     N_spike2 = T * rate2 |> Int

#     spikes1 = rand(N_spike1) * T
#     spikes2 = rand(N_spike2) * T
#     for n in eachindex(spikes1)
#         N = rand(1:length(spikes2))
#         spikes2[N] = spikes1[n] + τ_cov * randn()
#     end
#     N1 = fill([1], length(spikes1))
#     N2 = fill([2], length(spikes2))
#     neurons = vcat(N1, N2)
#     spiketimes = vcat(spikes1, spikes2)
#     inputs = SpikeTimeParameter(spiketimes, neurons)
#     st = Identity(N = max_neurons(inputs))
#     stim = SpikeTimeStimulusIdentity(st, :g, param = inputs)
#     w = zeros(Float32, 2, 2)
#     w[2, 1] = 5.0f0
#     syn = SpikingSynapse(st, st, :h, w = w, param = stdp_param)
#     model = compose(st = st, stim = stim, syn = syn, silent = true)
#     SNN.monitor!(model.pop..., [:fire])
#     train!(model = model, duration = T, dt = 0.1ms)
#     return model
#     ΔW = model.syn[1].W[1] - 5.0f0
#     return ΔW
# end

# export stdp_kernel, stdp_integral, stdp_weight_correlated

# function stdp_pairing()
#     plot(ylims = (-4, 4), xlims = (-10, 10), legend = false, frame = :none)
#     hline!([-2, 2], lc = :black, lw = 4)
#     scatter!([-6], [2.4], markershape = :vline, mc = :black, ms = 18, lw = 4)
#     scatter!([1], [-1.6], markershape = :vline, mc = :black, ms = 18, lw = 4)
#     annotate!(6, 2.4, text("Pre synaptic", 15, :bottom))
#     annotate!(6, -1.6, text("Post synaptic", 15, :bottom))
#     annotate!(-2.5, 0, text("ΔT", 15, :bottom))
#     plot!([-6, 1], [0, 0], lc = :black, arrow = (:both, :closed, 3))
# end

# ## Measure the weight change for decorrelated spike trains
# function stdp_weight_decorrelated(stdp_param, rate1 = 10Hz, rate2 = 10Hz)
#     st1 = Poisson(N = 50, param = PoissonParameter(rate = rate1))
#     st2 = Poisson(N = 50, param = PoissonParameter(rate = rate2))
#     syn = SpikingSynapse(st1, st2, nothing, p = 1.0f0, μ = 20, LTPParam = stdp_param)
#     model = compose(; st1, st2, syn, silent = true)
#     T = 20_000ms
#     train!(model = model, duration = T, dt = 0.1ms)
#     return (model.syn.syn.W .- 20) / T * 60^2
# end

# #


export stdp_kernel, stdp_kernel!
