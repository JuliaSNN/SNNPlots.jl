module SNNPlots

using ColorSchemes
using LaTeXStrings
using Measures
using SNNModels
import SNNModels: AbstractPopulation, AbstractStimulus, AbstractConnection
using UnPack
using Parameters
using Requires
@load_units

# @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("backend/makie.jl")
# @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" 
include("backend/plots.jl")
include(joinpath(@__DIR__,"plot.jl"))
include(joinpath(@__DIR__,"extra_plots.jl"))
include(joinpath(@__DIR__,"stdp_plots.jl"))
include(joinpath(@__DIR__,"spatial.jl"))
# end

export raster,
    vecplot,
    plot,
    plot!,
    save_model,
    load_model,
    plot_model,
    plot_stimulus,
    plot_connections
end
