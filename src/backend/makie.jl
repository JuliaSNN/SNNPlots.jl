using Makie


# onecolumn = (8.8cm, 13cm)
# twocolumn = (18cm, 185cm)

_backend = :Makie
@info "Using Makie backend for plotting"

const okabe_ito_10 = ColorScheme(get(getfield(ColorSchemes, :okabe_ito), range(0.0, 1.0, length=10)))[:]

macro makie_default()
    default_values = (;
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

    ex = Expr(:block, :(SNNPlots.set_theme!(;$(default_values)...)))
    esc(ex)
end

@makie_default
export @makie_default, default_colors, inch, cm, pt, nature_figure, okabe_ito_10
