using Makie


# onecolumn = (8.8cm, 13cm)
# twocolumn = (18cm, 185cm)

_backend = :Makie
@info "Using Makie backend for plotting"

"""
    okabe_ito_10

Ten colours sampled uniformly from the Okabe-Ito colour-blind-safe scheme
(`ColorSchemes.okabe_ito`), as a `Vector` of RGB colours. It is the default Makie colour cycle set
by [`makie_default!`](@ref) and is used by the spatial plots of SNNPlots.

# Example
```julia
using SpikingNeuralNetworks
c = okabe_ito_10[2]
```
"""
const okabe_ito_10 = ColorScheme(get(getfield(ColorSchemes, :okabe_ito), range(0.0, 1.0, length=10)))[:]

"""
    makie_default!()

Apply the SNNPlots Makie theme to the global Makie theme with `Makie.set_theme!`.

The theme sets: figure size `(600, 400)`; `Axis` without x and y grid lines; `Legend` without
frame and with label size 10; colour palette [`okabe_ito_10`](@ref); `Lines` and `Band` cycling
through the palette colours.

Since SNNPlots 0.2.10 it is called from `SNNPlots.__init__`, so it takes effect every time the
package is loaded (a top-level call would run only during precompilation and be lost in the user
session). Call it again to restore the theme after `set_theme!`/`update_theme!` elsewhere.
Returns the value of `set_theme!` (`nothing`).

# Example
```julia
using SpikingNeuralNetworks
SNNPlots.makie_default!()
```
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
"""
    @makie_default

Backward-compatible macro form of [`makie_default!`](@ref): expands to
`SNNPlots.makie_default!()`, i.e. re-applies the SNNPlots Makie theme. The theme is already applied
when SNNPlots is loaded, so the macro is only needed to restore it after another `set_theme!`.
The expansion refers to the name `SNNPlots`, which must therefore be visible in the calling scope
(it is, after `using SpikingNeuralNetworks` or `import SNNPlots`).
"""
macro makie_default()
    esc(:(SNNPlots.makie_default!()))
end

# `inch` and `pt` are the `Measures` lengths; `cm` resolves to the SNNModels length unit
# (`1cm == 1.0f0`, defined by `@load_units` in SNNPlots.jl), not to `Measures.cm`.
export @makie_default, makie_default!, inch, cm, pt, okabe_ito_10
