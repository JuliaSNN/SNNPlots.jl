# SNNPlots changelog

## 2.0.0

Requires SNNModels 2 (compat `SNNModels = "2"`); version number follows SNNModels.

- STDP plots read plastic weights through `plasticity_value` and build the synapse transpose arrays
  (`ensure_transpose!`) before reading `J`/`index`.
- Removed code that no module included: the Plots.jl backend (`backend/plots.jl`), `extra_plots.jl`
  (used an undefined `SNN.plasticity_value`), `old_plot.jl`, `other_plots.jl`. They were not loaded in 0.2.11
  either; the Makie backend is the only backend.
- Dependencies removed: Plots, DrWatson, Colors, Requires, LaTeXStrings, Parameters; CairoMakie is a test-only
  dependency (load a Makie backend, for example `using CairoMakie`, before plotting).
