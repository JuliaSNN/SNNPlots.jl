using Plots

_backend = :Plots
@info "Using Plots backend for plotting"

macro plots_default()
    default_values = (;
        bg_color_legend = :transparent,
        ytickfontsize = 12,
        xtickfontsize = 12,
        legend_title_font_halign = :right,
        legend_title_font_pointsize = 14,
        legend_font_pointsize = 11,
        margins = 5Plots.mm,
        palette = :okabe_ito,
        size = (900, 600),
        frame=:box, 
        grid=false, 
        tickfontsize=12, 
        guidefontsize=14, 
        legendfontsize=12, 
        margin=5Plots.mm, 
        foreground_color_legend=:transparent, 
        mscolor=:auto, 
        lw=4
    )

    ex = Expr(:block, :(default(;$(default_values)...)))
    esc(ex)
end

@plots_default
export plots_default