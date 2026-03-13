include("../tripod_network.jl")
# include("../experiments/analysis/spike_analysis.jl")
net = NetParams(
    neurons = 2500,
    types = duarte_types,
    model = "random",
    exc_noise = 1.5f0,
    sst_noise = 0.5f0,
    pv_noise = 0.5f0,
    conn = lkd2014,
)

stim = StimParams(
    density = 5.0f-2,
    simtime = 60000,
    duration = 50,
    input = "asymmetric",
    strength = 1.78f0,
    rate = 4.0f0,
    dictionary = "simple",
    mask_words = false,
)
store = StoreParams(id = joinpath("test", "plot"), interval = 1000)
println("Data stored in ", store.data)


seq = seq_encoder(net, stim)
learn = LearnParams()
dends = get_dends(net)
store_network_params(stim, seq, net, dends, learn, store)

##

tripod_network(
    store.path,
    store_weights = true,
    store_spikes = true,
    store_membrane = true,
    debug = true,
)
##
using ColorSchemes
params = load_params(store.path)
sequence = params.seq.sequence
word = params.seq.populations[sequence[1, 6]]
cc = palette(:roma, 5)
q = plot(save_v[1, word, 1:1:7500]', labels = "", c = :black, ylims = (-100, 20))
word = params.seq.mapping[sequence[1, 6]]
annotate!((2800, 15, Plots.text(word, :black, 11)))
vline!([500 * 4], ls = :dash, c = :red, label = "")
vline!([500 * 7], ls = :dash, c = :red, label = "")

save_v


plots = []
for x = 5:7
    seq_n = sequence[2, x]
    ph = params.seq.populations[seq_n]
    active1 = findall(params.seq.dendrites[seq_n][1, :] .> 0)
    active2 = findall(params.seq.dendrites[seq_n][2, :] .> 0)
    p = plot(
        save_v[2, ph[active1], 1:1:7500]',
        labels = "",
        ylims = (-100, 10),
        c = cc[x-2],
        alpha = 0.5,
    )
    p = plot!(
        save_v[3, ph[active2], 1:1:7500]',
        labels = "",
        ylims = (-100, 10),
        c = cc[x-2],
        alpha = 0.5,
    )
    # p = plot()
    plot!(
        p,
        (
            mean(save_v[2, ph, 1:1:7500], dims = 1)[1, :] .+
            mean(save_v[3, ph, 1:7500], dims = 1)[1, :]
        ) ./ 2,
        c = :black,
        labels = "",
        xticks = :none,
        xaxis = false,
    )
    push!(plots, p)
    ph = params.seq.mapping[sequence[2, x]]
    vline!([500 * (x - 1)], ls = :dash, c = :black, label = "")
    vline!([500 * (x)], ls = :dash, c = :black, label = "")
    annotate!((100 + 500 * (x - 1), 1, Plots.text(ph, :black, 11, :left)))
end
# plot!(plots[3], ylabel="      Membrane potential")
plot!(q, xlabel = "Time (ms)")
pp = plot!(plots..., q, layout = (4, 1), yticks = :none)


##
inh_spikes = 0.0f0
exc_spikes = 0.0f0
X = 100.0f0

g = zeros(Float32, 4, 4, 2)
t = 2500
rec = zeros(Float32, t, 4, 4, 2)
cc = 1
for tt = 1:t
    if (tt == 100) || tt == 300 || tt == 500
        inh_spikes = X
        exc_spikes = X
    else
        inh_spikes = 0.0f0
        exc_spikes = 0.0f0
    end
    @views update_synapse_soma!(g[1, :, :], inh_spikes, exc_spikes, Esyn_soma)
    @views update_synapse_dend!(g[2, :, :], inh_spikes, exc_spikes, Esyn_dend)
    @views update_synapse_dend!(g[3, :, :], inh_spikes, exc_spikes, Isyn_sst)
    @views update_synapse_dend!(g[4, :, :], inh_spikes, exc_spikes, Isyn_pv)
    rec[tt, :, :, :] .= g[:, :, :]
end

rec[:, 1, 1, 1]
##
##
pcond = begin
    p = plot()
    ll = ["AMPA", "NMDA", "GABAa", "GABAb"]
    for (n, name) in enumerate(fieldnames(Receptors))
        try
            g = [getfield(getfield(Esyn_soma, name), :gsyn)]
            bar!([n - 0.3], g, c = :orange, barwidths = 0.2, label = "")
            g = [getfield(getfield(Esyn_dend, name), :gsyn)]
            bar!([n - 0.1], g, c = :green, barwidths = 0.2, label = "")
            g = [getfield(getfield(Isyn_sst, name), :gsyn)]
            bar!([n + 0.1], g, c = :darkred, barwidths = 0.2, label = "")
            g = [getfield(getfield(Isyn_pv, name), :gsyn)]
            bar!([n + 0.3], g, c = :lightblue, barwidths = 0.2, label = "")
            plot!(xticks = (1:4, ll), ylabel = "Conductance", guidefontsize = 18)
        catch
            scatter!(
                [[], [], [], []],
                c = [:orange :green :darkred :lightblue],
                label = ["Soma" "Dend" "SST" "PV"],
                msc = [:orange :green :darkred :lightblue],
                legendfontsize = 12,
            )
            continue
        end
    end
    p
end


tt = ["soma", "dend", "sst", "pv"]
ll = [["AMPA" "NMDA" "GABAa" "GABAb"], "", "", ""]
psynapse = plot(
    [plot(rec[90:end, x, :, 2], title = tt[x], label = ll[x]) for x = 1:4]...,
    guidefontsize = 18,
    xticks = (0:500:2500, 0:50:250),
    legendfontsize = 10,
    titlefontsize = 14,
    tickfontsize = 8,
)

layout = @layout [
    a{0.4h}
    b{0.6h}
]
savefig(
    plot(pcond, psynapse, layout = layout),
    joinpath(@__DIR__, "receptor_conductances.pdf"),
)
