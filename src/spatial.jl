"""
    plot_spatial_connectivity(fig, connectivity, config; ...)

Visualizes the spatial layout of neuronal populations and their connections to a single postsynaptic neuron.

Plots all neurons from Exc, PV, SST populations with distinct markers. For a specified postsynaptic
neuron (`n0` from population `post`), highlights all presynaptic neurons connected to it and draws
arrows indicating those connections.

# Arguments
- `fig`: Makie figure or GridLayout position to plot into.
- `connectivity`: structure containing neuron `points` (locations) and `links` (adjacency matrices).
- `config`: network configuration, containing `network` and `spatial` properties.
- `post::Symbol`: postsynaptic population (default: `:Exc`).
- `n0::Int`: index of the target postsynaptic neuron.
- `do_arrows::Bool`: draw connection arrows (default: `true`).
- `rotation::Symbol`: `:horizontal` (default, x = tonotopic) or `:vertical` (swap x↔y axes).
- `xlabel`, `ylabel`: axis labels for the x-data and y-data axes; default `""`. Swapped automatically
  when `rotation = :vertical`. The Legend is unaffected by rotation.
"""
function plot_spatial_connectivity(fig, connectivity, config;
        post      = :Exc,
        n0::Int   = rand(1:20),
        do_arrows = true,
        rotation::Symbol = :horizontal,
        xlabel    = "",
        ylabel    = "",
        do_legend = true)

    @unpack points, links = connectivity
    @unpack network, spatial = config
    @unpack recurrence, Npop = network
    @unpack connections = recurrence
    @unpack grid_size = spatial

    colors = okabe_ito_10[[8, 2, 3]]
    shapes = [:utriangle, :circle, :hexagon]
    ticks  = (range(0, 0.1, 5), string.(range(0, 1, 5)))

    # swap x↔y for :vertical; identity for :horizontal
    pt = rotation == :vertical ? (x, y) -> Point2f(y, x) : (x, y) -> Point2f(x, y)
    ax_xlabel = rotation == :vertical ? ylabel : xlabel
    ax_ylabel = rotation == :vertical ? xlabel : ylabel

    xyns, cs, ms_list = Point2f[], Any[], Symbol[]
    ax = Axis(fig[1:2, 1:2];
            #   aspect  = DataAspect(),
              xlabel  = ax_xlabel,
              ylabel  = ax_ylabel,
              xticks  = ticks,
              yticks  = ticks)

    for (pre, c, ms) in zip([:Exc, :PV, :SST], colors, shapes)
        w = []
        for n in eachindex(points[pre])
            conn    = network.recurrence.connections[name(pre, post)]
            targets = haskey(conn, :target) ? config.network.targets[conn.target] : [nothing]
            comp    = targets[1]
            _name   = str_name(pre, post, comp)
            push!(w, links[Symbol(_name)][n0, n])
        end

        x  = [points[pre][i][1] for i in eachindex(points[pre])]
        y  = [points[pre][i][2] for i in eachindex(points[pre])]
        xy = pt.(x, y)
        xyn = pt.(x[findall(v -> v == 1, w)], y[findall(v -> v == 1, w)])
        append!(xyns, xyn)
        append!(cs, fill(c, length(xyn)))
        append!(ms_list, fill(ms, length(xyn)))

        Makie.scatter!(ax, xy; color = c, markersize = 2, alpha = 0.5, marker = ms)
    end

    xy0 = pt(points[post][n0][1], points[post][n0][2])

    for n in eachindex(xyns)
        if do_arrows
            v = xy0 .- xyns[n]
            Makie.arrows2d!(ax, xyns[n], v; shaftwidth = 1, tiplength = 0, color = cs[n], alpha = 0.1)
        end
        Makie.scatter!(ax, xyns[n]; color = cs[n], markersize = 5, marker = ms_list[n], strokewidth = 0.1)
    end
    Makie.scatter!(ax, xy0; marker = :cross, markersize = 15, strokecolor = :black, strokewidth = 2, color = :white)

    markers = map(zip(shapes, colors)) do (ms, c)
        MarkerElement(color = c, marker = ms, markersize = 10)
    end
    if do_legend 
        Legend(fig[0,1:2], markers, ["Exc", "PV", "SST"],
            ""; position = :rt, orientation = :horizontal,
            tellheight = false, tellwidth = false)
        return ax
    else
    legend_info = (; markers, labels = ["Exc", "PV", "SST"], position = :rt, orientation = :horizontal, tellheight = false, tellwidth = false)
    return (ax, legend_info)
    end
end


"""
    plot_connection_distances(fig; ds, rs)

Computes and plots the distribution of connection distances for different presynaptic populations
to the excitatory (`:Exc`) population.

This function calculates the periodic distance for every established connection from presynaptic
populations (`:Exc`, `:PV`, `:SST`) to all postsynaptic `:Exc` neurons. It then generates two plots:
1. A bar plot showing the binned connection counts (density).
2. A line plot showing the connection probability, normalized by the area of the annulus at each distance.

# Arguments
- `fig`: The Makie `Figure` object to plot into.
- `ds`: A dictionary containing the binned connection counts for each presynaptic population.
- `rs`: The edges of the distance bins used for histogramming the connection distances.
"""
function plot_connection_distances(fig; ds, rs, probability = true)

    ax1 = Axis(fig[1,1], xlabel="Distance (mm)", ylabel="Conn. density", )
    if probability
        ax2 = Axis(fig[1,2], xlabel="Distance (mm)", ylabel="Conn. probability")
    end
    post = :Exc
    scaled_ws = Float32[]
    ws = Float32[]
    grp = Int[]
    # dodge = []
    xs = rs[2:end].-Float32(rs.step/2)
    colors = okabe_ito_10[[8, 2, 3]]
    for (n, pre) in enumerate([:Exc, :PV, :SST])
        w = ds[name(pre, post)] 
        append!(ws, copy(w./sum(w)))
        w = w./rs[2:end].^2
        if probability
            lines!(ax2, xs.*10, w./sum(w) , label = string(pre), color = colors[n], linewidth=4)
        end
        append!(scaled_ws, copy(w./sum(w)))
        append!(grp, fill(n, length(w)))
    end

    xxs = repeat(xs, outer=3)
    grp
    colors = repeat(okabe_ito_10[[8, 2, 3]], inner=length(xs))
    barplot!(ax1, xxs.*10, ws, dodge=grp, color = colors, label = ["Exc", "PV", "SST"])

    colors = okabe_ito_10[[8, 2, 3]]
    # markers = map(zip(shapes, colors)) do (ms, c)
    #     PolyElement(;color=c)
    # end
    # axislegend(ax2, position = :rt)
    # axislegend(ax1, markers, ["Exc", "PV", "SST"], position = :rt)
end

export plot_spatial_connectivity, plot_connection_distances
#

# function plot_connected_neurons(
#     connectivity,
#     config,
#     post::Symbol = :Exc,
#     n0::Int = rand(1:20),
# )
#     @unpack points, links = connectivity
#     @unpack network, spatial = config
#     @unpack connections, Npop = network
#     @unpack grid_size = spatial
#     plots = map(zip([:Exc, :PV, :SST], [:darkred, :blue, :darkorange])) do (pre, c)
#         d = []
#         w = []
#         for n = 1:length(points[pre])
#             # push!(d, periodic_distance(points[post][n0], points[pre][n], grid_size))
#             conn = network.connections[name(pre, post)]
#             targets = haskey(conn, :target) ? conn.target : [nothing]
#             comp = targets[1]
#             _name = str_name(pre, post, comp)
#             ll = links[Symbol(_name)][n0, n]
#             push!(w, ll)
#             # ll && continue 
#         end

#         xlims = extrema([
#             getfield(points, post)[i][1] for i in eachindex(getfield(points, post))
#         ])
#         ylims = extrema([
#             getfield(points, post)[i][2] for i in eachindex(getfield(points, post))
#         ])
#         x = [points[pre][i][1] for i = 1:length(points[pre])]
#         y = [points[pre][i][2] for i = 1:length(points[pre])]
#         x0 = points[post][n0][1]
#         y0 = points[post][n0][2]
#         scatter(
#             x,
#             y,
#             xlabel = "X",
#             ylabel = "Y",
#             title = "",
#             c = :grey,
#             ms = 2,
#             msc = :grey,
#             label = "",
#             alpha = 0.4,
#         )
#         n = findall(x->x==1, w)
#         scatter!(x[n], y[n], c = c, ms = 5, msc = c, label = "")
#         scatter!(
#             [x0],
#             [y0],
#             xlabel = "Tonotopic axis (mm)",
#             ylabel = "Isofrequency axis (mm)",
#             title = "Connected neurons",
#             c = :black,
#             ms = 8,
#             msc = :grey,
#             label = "",
#         )
#         ticks = (range(0, 0.1, 5), range(0, 1, 5))
#         plot!(
#             title = string("$pre => $post"),
#             xlims = xlims,
#             ylims = ylims,
#             frame = :axes,
#             margin = 10Plots.mm,
#             xticks = ticks,
#             yticks = ticks,
#         )
#     end
#     p1 = plot(plots..., size = (1200, 400), layout = (1, 3))


#     plots = map([:Exc, :PV, :SST]) do pre
#         conn = network.connections[name(pre, post)]
#         targets = haskey(conn, :target) ? conn.target : [nothing]
#         comp = targets[1]
#         _name = str_name(pre, post, comp) |> Symbol
#         _m = maximum([maximum(sum(links[_name], dims = 2)[:, 1]), 3])
#         p = create_connection_histogram(connectivity, config, pre, post)
#         # N_s, N_l = count_neurons(pre, post, network, config)
#         # annotate!(p,(0.1,0.8), text("Short: $(N_s)/$(N_s+N_l) \nLong: $(N_l)/$(N_s+N_l)", 12, :left, :white))
#     end
#     p2 = plot(plots..., size = (1200, 400), layout = (1, 3))
#     plot(p1, p2, layout = (2, 1), size = (1200, 800))
# end

# function count_neurons(pre::Symbol, post::Symbol, network, config, samples = 50)
#     @unpack connections, Npop, spatial = config
#     @unpack grid_size, ϵ, p_long = spatial
#     @unpack points, links = network
#     _Ns = 0
#     _Nl = 0
#     for n0 in rand(1:Npop.Exc, samples)
#         for n in eachindex(points[pre])
#             !(links[name(pre, post)][n0, n]) && continue
#             _d = periodic_distance(points[post][n0], points[pre][n], grid_size)
#             if _d < 0.2
#                 _Ns += 1
#             else
#                 _Nl += 1
#             end
#         end
#     end
#     _Ns, _Nl = round(Int, _Ns/samples), round(Int, _Nl/samples)
#     N_s =
#         ϵ * connections[name(pre, post)].p * Npop[pre] * (1-p_long[pre]) |>
#         x -> round(Int, x)
#     N_l =
#         ϵ * connections[name(pre, post)].p * Npop[pre] * (p_long[pre]) |> x -> round(Int, x)
#     @info "$pre => $post (real/expected):  Short ($_Ns/$(N_s));  Long ($(_Nl)/$(N_l))"
#     return N_s, N_l
# end

# # Function to create a 2D histogram of connection probabilities
# function create_connection_histogram(connectivity, config, pre = :Exc, post = :Exc)

#     @unpack points, links = connectivity
#     @unpack grid_size = config.spatial
#     @unpack connections = config

#     grid_size = typeof(grid_size) <: Real ? [grid_size, grid_size] : grid_size
#     @assert length(grid_size) <= 2 "Grid size must be a 2D vector [width, height]"

#     # Initialize the histogram
#     bins = 100
#     hist = zeros(bins, bins)
#     counts = zeros(bins, bins)

#     pre_points = getfield(points, pre)
#     post_points = getfield(points, post)
#     xlims = extrema([points[i][1] for i in eachindex(points)])
#     ylims = extrema([points[i][1] for i in eachindex(points)])

#     for j in eachindex(pre_points)
#         conn = config.network.connections[name(pre, post)]
#         targets = haskey(conn, :target) ? conn.target : [nothing]
#         comp = targets[1]
#         _name = str_name(pre, post, comp) |> Symbol
#         for i in findall(links[_name][:, j])
#             # distance = periodic_distance(post_points[i], pre_points[j], grid_size)
#             x = post_points[i][1] - pre_points[j][1]
#             y = post_points[i][2] - pre_points[j][2]
#             x = min(abs(x), grid_size[1] - abs(x))*sign(x)
#             y = min(abs(y), grid_size[2] - abs(y))*sign(y)
#             bin_x = Int(floor(x * bins / grid_size[1])) + bins ÷ 2 + 1
#             bin_y = Int(floor(y * bins / grid_size[2])) + bins ÷ 2 + 1
#             counts[bin_x, bin_y] += 1
#         end
#     end
#     heatmap(
#         range(0, 1, bins),
#         range(0, 1, bins),
#         counts',
#         xlabel = "Distance",
#         ylabel = "Distance",
#         title = "",
#         color = :viridis,
#         cbar = false,
#     )
# end


# # Function to create the animation
# function animate_raster(points, model, interval, path, fps = 10)
#     Exc_points = points.Exc
#     PV_points = points.PV
#     SST_points = points.SST

#     anim = @animate for t in interval
#         plot(
#             xlim = (0, 1),
#             ylim = (0, 1),
#             size = (800, 800),
#             legend = :none,
#             title = "t = $t ms",
#         )
#         # for (pop, points, c) in zip([:PV, :SST, :Exc], [PV_points, SST_points, Exc_points], [:blue, :orange, :black])
#         begin
#             for delay = 1:2:10
#                 interval = [t-delay, t-delay+2SNN.ms]
#                 c = :green
#                 _, n = SNN._raster(spiketimes(model.pop[:PV]), interval)
#                 n = Int.(n)
#                 scatter!(
#                     [p[1] for p in PV_points[n]],
#                     [p[2] for p in PV_points[n]],
#                     color = c,
#                     msc = c,
#                     markersize = 15,
#                     alpha = 0.3/log(delay+1),
#                 )
#                 c = :orange
#                 _, n = SNN._raster(spiketimes(model.pop[:SST]), interval)
#                 n = Int.(n)
#                 scatter!(
#                     [p[1] for p in SST_points[n]],
#                     [p[2] for p in SST_points[n]],
#                     color = c,
#                     msc = c,
#                     markersize = 15,
#                     alpha = 0.3/log(delay+1),
#                 )
#             end
#         end
#         pop = :Exc
#         points = Exc_points
#         c = :black
#         begin
#             for delay = 1:2:30
#                 interval = [t-delay, t-delay+2ms]
#                 _, n = SNN._raster(spiketimes(model.pop[pop]), interval)
#                 n = Int.(n)
#                 scatter!(
#                     [p[1] for p in points[n]],
#                     [p[2] for p in points[n]],
#                     color = c,
#                     msc = c,
#                     markersize = 5/(delay+0.1),
#                     alpha = 1/(delay+0.1),
#                 )
#             end

#             _, n = SNN._raster(spiketimes(model.pop[pop]), [t, t+2ms])
#             n = Int.(n)
#             scatter!(
#                 [p[1] for p in points[n]],
#                 [p[2] for p in points[n]],
#                 color = c,
#                 msc = c,
#                 markersize = 6,
#             )
#         end
#     end

#     gif(anim, path, fps = 10)
# end

# export plot_connected_neurons, create_connection_histogram, animate_raster
