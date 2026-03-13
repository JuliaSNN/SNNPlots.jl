using Makie


# onecolumn = (8.8cm, 13cm)
# twocolumn = (18cm, 185cm)

function nature_figure(column::Int, ratio)
    inch = 96
    pt =4/3
    cm = inch/2.54
    if column == 1
        width, height = 8.8cm, 22.0cm
    elseif column == 2
        width, height = 18cm, 22.5cm
    else
        error("Column must be 1 or 2")
    end
    return (width, height*ratio)
end

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
