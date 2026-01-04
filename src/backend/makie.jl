using Makie

_backend = :Makie
@info "Using Makie backend for plotting"

set_theme!(
    # theme_dark(),
    fontsize = 20,
    palette = (color = 
    get(ColorSchemes.viridis, range(0.5, 1.0, length=2)) |> ColorScheme
    , linestyle = [:dash, :dot])
)