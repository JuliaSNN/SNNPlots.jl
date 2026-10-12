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
    using Measures
    using SNNModels
    import SNNModels: AbstractPopulation, AbstractStimulus, AbstractConnection
    using UnPack

    _backend = nothing
    function __init__()
        # the theme must be set at load time, not at precompile time (see makie_default!)
        makie_default!()
    end

    include("backend/makie.jl")
    include(joinpath(@__DIR__,"raster.jl"))
    include(joinpath(@__DIR__,"vecplot.jl"))
    include(joinpath(@__DIR__,"stdp_plots.jl"))
    include(joinpath(@__DIR__,"spatial.jl"))

    @load_units
    # `plot`, `plot!` are re-exported from Makie and `save_model`, `load_model` from SNNModels.
    export raster,
        raster!,
        vecplot,
        plot,
        plot!,
        save_model,
        load_model,
        stdp_kernel,
        stdp_kernel!

end
