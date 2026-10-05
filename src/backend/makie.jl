using Makie


# onecolumn = (8.8cm, 13cm)
# twocolumn = (18cm, 185cm)

_backend = :Makie
@info "Using Makie backend for plotting"

const okabe_ito_10 = ColorScheme(get(getfield(ColorSchemes, :okabe_ito), range(0.0, 1.0, length=10)))[:]

"""
    makie_default!()

Apply the SNNPlots Makie theme (no grid, frameless legends, Okabe-Ito palette) to the global
Makie theme. Called from `SNNPlots.__init__`, so it takes effect every time the package is
loaded; a top-level call would run only during precompilation and be lost in the user session.
Call it again to restore the theme after `set_theme!`/`update_theme!` elsewhere.
"""
function makie_default!()
    set_theme!(;
        size = (600, 400),
        # fontsize = 10pt,
        # yticklabelsize = 8pt,
        # xticklabelsize = 8pt,
        Axis = (
            xgridvisible = false,
            ygridvisible = false,
        ),
        Legend = (;
            framevisible = false,
            labelsize = 10,
        ),
        palette = (;color =  okabe_ito_10
        ),
        Band = (;cycle =:color),
        Lines = (;cycle =:color),
    )
end

# kept for backward compatibility: `@makie_default` == `SNNPlots.makie_default!()`
macro makie_default()
    esc(:(SNNPlots.makie_default!()))
end

export @makie_default, makie_default!, default_colors, inch, cm, pt, nature_figure, okabe_ito_10
