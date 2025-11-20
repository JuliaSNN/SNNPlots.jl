module SNNPlots

using Plots
using ColorSchemes
using LaTeXStrings
using Measures
using SNNModels
import SNNModels: AbstractPopulation, AbstractStimulus, AbstractConnection
using UnPack
using Parameters
@load_units

include("plot.jl")
include("extra_plots.jl")
include("stdp_plots.jl")
include("spatial.jl")

# default(fg_legend = :transparent)
default(bg_color_legend = :transparent)
default(xguidefontsize = 15)
default(yguidefontsize = 15)
default(ytickfontsize = 12)
default(xtickfontsize = 12)
default(legend_title_font_halign = :right)
default(legend_title_font_pointsize = 14)
default(legend_font_pointsize = 11)
default(margins = 5Plots.mm)


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
