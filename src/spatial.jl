
function plot_connected_neurons(
    connectivity,
    config,
    post::Symbol = :Exc,
    n0::Int = rand(1:20),
)
    @unpack points, links = connectivity
    @unpack network, spatial = config
    @unpack connections, Npop = network
    @unpack grid_size = spatial
    plots = map(zip([:Exc, :PV, :SST], [:darkred, :blue, :darkorange])) do (pre, c)
        d = []
        w = []
        for n = 1:length(points[pre])
            # push!(d, periodic_distance(points[post][n0], points[pre][n], grid_size))
            conn = network.connections[name(pre, post)]
            targets = haskey(conn, :target) ? conn.target : [nothing]
            comp = targets[1]
            _name = str_name(pre, post, comp)
            ll = links[Symbol(_name)][n0, n]
            push!(w, ll)
            # ll && continue 
        end

        xlims = extrema([
            getfield(points, post)[i][1] for i in eachindex(getfield(points, post))
        ])
        ylims = extrema([
            getfield(points, post)[i][2] for i in eachindex(getfield(points, post))
        ])
        x = [points[pre][i][1] for i = 1:length(points[pre])]
        y = [points[pre][i][2] for i = 1:length(points[pre])]
        x0 = points[post][n0][1]
        y0 = points[post][n0][2]
        scatter(
            x,
            y,
            xlabel = "X",
            ylabel = "Y",
            title = "",
            c = :grey,
            ms = 2,
            msc = :grey,
            label = "",
            alpha = 0.4,
        )
        n = findall(x->x==1, w)
        scatter!(x[n], y[n], c = c, ms = 5, msc = c, label = "")
        scatter!(
            [x0],
            [y0],
            xlabel = "Tonotopic axis (mm)",
            ylabel = "Isofrequency axis (mm)",
            title = "Connected neurons",
            c = :black,
            ms = 8,
            msc = :grey,
            label = "",
        )
        ticks = (range(0, 0.1, 5), range(0, 1, 5))
        plot!(
            title = string("$pre => $post"),
            xlims = xlims,
            ylims = ylims,
            frame = :axes,
            margin = 10Plots.mm,
            xticks = ticks,
            yticks = ticks,
        )
    end
    p1 = plot(plots..., size = (1200, 400), layout = (1, 3))


    plots = map([:Exc, :PV, :SST]) do pre
        conn = network.connections[name(pre, post)]
        targets = haskey(conn, :target) ? conn.target : [nothing]
        comp = targets[1]
        _name = str_name(pre, post, comp) |> Symbol
        _m = maximum([maximum(sum(links[_name], dims = 2)[:, 1]), 3])
        p = create_connection_histogram(
            connectivity,
            config.spatial,
            config.network,
            pre,
            post,
        )
        # N_s, N_l = count_neurons(pre, post, network, config)
        # annotate!(p,(0.1,0.8), text("Short: $(N_s)/$(N_s+N_l) \nLong: $(N_l)/$(N_s+N_l)", 12, :left, :white))
    end
    p2 = plot(plots..., size = (1200, 400), layout = (1, 3))
    plot(p1, p2, layout = (2, 1), size = (1200, 800))
end

function count_neurons(pre::Symbol, post::Symbol, network, config, samples = 50)
    @unpack connections, Npop, spatial = config
    @unpack grid_size, ϵ, p_long = spatial
    @unpack points, links = network
    _Ns = 0
    _Nl = 0
    for n0 in rand(1:Npop.Exc, samples)
        for n in eachindex(points[pre])
            !(links[name(pre, post)][n0, n]) && continue
            _d = periodic_distance(points[post][n0], points[pre][n], grid_size)
            if _d < 0.2
                _Ns += 1
            else
                _Nl += 1
            end
        end
    end
    _Ns, _Nl = round(Int, _Ns/samples), round(Int, _Nl/samples)
    N_s =
        ϵ * connections[name(pre, post)].p * Npop[pre] * (1-p_long[pre]) |>
        x -> round(Int, x)
    N_l =
        ϵ * connections[name(pre, post)].p * Npop[pre] * (p_long[pre]) |> x -> round(Int, x)
    @info "$pre => $post (real/expected):  Short ($_Ns/$(N_s));  Long ($(_Nl)/$(N_l))"
    return N_s, N_l
end

# Function to create a 2D histogram of connection probabilities
function create_connection_histogram(
    connectivity,
    spatial,
    network,
    pre = :Exc,
    post = :Exc,
)

    @unpack points, links = connectivity
    @unpack grid_size = spatial
    @unpack connections = network

    grid_size = typeof(grid_size) <: Real ? [grid_size, grid_size] : grid_size
    @assert length(grid_size) <= 2 "Grid size must be a 2D vector [width, height]"

    # Initialize the histogram
    bins = 100
    hist = zeros(bins, bins)
    counts = zeros(bins, bins)

    pre_points = getfield(points, pre)
    post_points = getfield(points, post)
    xlims = extrema([points[i][1] for i in eachindex(points)])
    ylims = extrema([points[i][1] for i in eachindex(points)])

    for j in eachindex(pre_points)
        conn = config.network.connections[name(pre, post)]
        targets = haskey(conn, :target) ? conn.target : [nothing]
        comp = targets[1]
        _name = str_name(pre, post, comp) |> Symbol
        for i in findall(links[_name][:, j])
            # distance = periodic_distance(post_points[i], pre_points[j], grid_size)
            x = post_points[i][1] - pre_points[j][1]
            y = post_points[i][2] - pre_points[j][2]
            x = min(abs(x), grid_size[1] - abs(x))*sign(x)
            y = min(abs(y), grid_size[2] - abs(y))*sign(y)
            bin_x = Int(floor(x * bins / grid_size[1])) + bins ÷ 2 + 1
            bin_y = Int(floor(y * bins / grid_size[2])) + bins ÷ 2 + 1
            counts[bin_x, bin_y] += 1
        end
    end
    heatmap(
        range(0, 1, bins),
        range(0, 1, bins),
        counts',
        xlabel = "Distance",
        ylabel = "Distance",
        title = "",
        color = :viridis,
        cbar = false,
    )
end


# Function to create the animation
function animate_raster(points, model, interval, path, fps = 10)
    Exc_points = points.Exc
    PV_points = points.PV
    SST_points = points.SST

    anim = @animate for t in interval
        plot(
            xlim = (0, 1),
            ylim = (0, 1),
            size = (800, 800),
            legend = :none,
            title = "t = $t ms",
        )
        # for (pop, points, c) in zip([:PV, :SST, :Exc], [PV_points, SST_points, Exc_points], [:blue, :orange, :black])
        begin
            for delay = 1:2:10
                interval = [t-delay, t-delay+2SNN.ms]
                c = :green
                _, n = SNN._raster(spiketimes(model.pop[:PV]), interval)
                n = Int.(n)
                scatter!(
                    [p[1] for p in PV_points[n]],
                    [p[2] for p in PV_points[n]],
                    color = c,
                    msc = c,
                    markersize = 15,
                    alpha = 0.3/log(delay+1),
                )
                c = :orange
                _, n = SNN._raster(spiketimes(model.pop[:SST]), interval)
                n = Int.(n)
                scatter!(
                    [p[1] for p in SST_points[n]],
                    [p[2] for p in SST_points[n]],
                    color = c,
                    msc = c,
                    markersize = 15,
                    alpha = 0.3/log(delay+1),
                )
            end
        end
        pop = :Exc
        points = Exc_points
        c = :black
        begin
            for delay = 1:2:30
                interval = [t-delay, t-delay+2ms]
                _, n = SNN._raster(spiketimes(model.pop[pop]), interval)
                n = Int.(n)
                scatter!(
                    [p[1] for p in points[n]],
                    [p[2] for p in points[n]],
                    color = c,
                    msc = c,
                    markersize = 5/(delay+0.1),
                    alpha = 1/(delay+0.1),
                )
            end

            _, n = SNN._raster(spiketimes(model.pop[pop]), [t, t+2ms])
            n = Int.(n)
            scatter!(
                [p[1] for p in points[n]],
                [p[2] for p in points[n]],
                color = c,
                msc = c,
                markersize = 6,
            )
        end
    end

    gif(anim, path, fps = 10)
end

export plot_connected_neurons, create_connection_histogram, animate_raster
