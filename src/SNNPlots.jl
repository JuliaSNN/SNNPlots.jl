"""
    SNNPlots

Makie-based plotting for SNNModels: raster plots ([`raster`](@ref), [`raster!`](@ref)),
recorded-variable traces ([`vecplot`](@ref), [`vecplot!`](@ref)), STDP kernels
([`stdp_kernel`](@ref), [`stdp_kernel!`](@ref), [`stdp_test`](@ref)) and spatial-network plots
([`plot_spatial_connectivity`](@ref), [`plot_connection_distances`](@ref)).

SNNPlots depends on `Makie` only; load a Makie backend (`using CairoMakie` or `using GLMakie`) to
display or save figures. The SNNPlots theme ([`makie_default!`](@ref)) is applied when the
package is loaded.
"""
module SNNPlots

    using ColorSchemes
    using LaTeXStrings
    using Measures
    using SNNModels
    import SNNModels: AbstractPopulation, AbstractStimulus, AbstractConnection
    using UnPack
    using Parameters
    using Requires

    _backend = nothing
    function __init__()
        # the theme must be set at load time, not at precompile time (see makie_default!)
        makie_default!()
        # @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("backend/plots.jl")
        # @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("backend/makie.jl")

        # if _backend == :Makie
        #     @info "Makie backend set with default settings"
        #     include(joinpath(@__DIR__,"raster.jl"))
        #     include(joinpath(@__DIR__,"vecplot.jl"))
        # end
        # if _backend == :Plots
        #     @info "Plots backend set with default settings"
        #     include(joinpath(@__DIR__,"raster.jl"))
        #     include(joinpath(@__DIR__,"vecplot.jl"))
        #     include(joinpath(@__DIR__,"other_plots.jl"))
        #     include(joinpath(@__DIR__,"extra_plots.jl"))
        #     include(joinpath(@__DIR__,"stdp_plots.jl"))
        #     include(joinpath(@__DIR__,"spatial.jl"))
        # end
    end

    # end
    include("backend/makie.jl")
    include(joinpath(@__DIR__,"raster.jl"))
    include(joinpath(@__DIR__,"vecplot.jl"))
    include(joinpath(@__DIR__,"stdp_plots.jl"))
    include(joinpath(@__DIR__,"spatial.jl"))

    @load_units
    # `plot`, `plot!` are re-exported from Makie and `save_model`, `load_model` from SNNModels;
    # `plot_model`, `plot_stimulus`, `plot_connections` are not defined in SNNPlots 0.2.10.
    export raster,
        vecplot,
        plot,
        plot!,
        save_model,
        load_model,
        plot_model,
        plot_stimulus,
        plot_connections,
        stdp_kernel,
        stdp_kernel!

end
