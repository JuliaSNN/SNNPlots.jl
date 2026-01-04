using Plots

_backend = :Plots
@info "Using Plots backend for plotting"

default(bg_color_legend = :transparent)
default(ytickfontsize = 12)
default(xtickfontsize = 12)
default(legend_title_font_halign = :right)
default(legend_title_font_pointsize = 14)
default(legend_font_pointsize = 11)
default(margins = 5Plots.mm)
default(palette = :okabe_ito)
default(size = (600, 400))
default(frame=:box,  
              grid=false, 
              tickfontsize=12, 
              guidefontsize=14, 
              legendfontsize=12, 
              margin=5Plots.mm, foreground_color_legend=:transparent, mscolor=:auto, lw=4)