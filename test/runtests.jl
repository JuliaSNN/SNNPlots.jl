using SNNPlots
using SNNModels
using CairoMakie
using Test
SNNModels.@load_units

@testset "SNNPlots" begin
    @testset "exports are defined" begin
        @test isempty([n for n in names(SNNPlots) if !isdefined(SNNPlots, n)])
        @test isdefined(SNNPlots, :raster!) && :raster! in names(SNNPlots)
    end

    @testset "stdp_test measures one pre/post pair" begin
        p = STDPGerstner()
        dw_plus = stdp_test(p; ΔT = 10ms)
        dw_minus = stdp_test(p; ΔT = -10ms)
        # Gerstner pair rule: ΔW(+Δ) = A_post' e^{-Δ/τ}, ΔW(-Δ) of opposite sign
        @test sign(dw_plus) == -sign(dw_minus)
        @test isapprox(abs(dw_plus), abs(dw_minus); rtol = 0.02)
        @test isapprox(abs(dw_plus), 1.0e-4 * exp(-10 / 20); rtol = 0.02)
    end

    E = SNNModels.Poisson(N = 20, param = PoissonParameter(20Hz), name = "E")
    I = SNNModels.Poisson(N = 10, param = PoissonParameter(30Hz), name = "I")
    monitor!([E, I], [:fire])
    sim!([E, I]; duration = 500ms)

    @testset "raster" begin
        fig, ax, plt = raster(spiketimes(E))
        @test ax.xlabel[] == "Time (ms)"
        fig, ax, plt = raster(spiketimes(E), [0, 500]; markersize = 3)
        f = Figure()
        ax1 = Axis(f[1, 1]); ax2 = Axis(f[2, 1])
        raster!(ax1, (; E, I), 0:500ms; names = ["exc", "inh"], marker = :rect)
        @test ax1.yticks[][2] == ["exc", "inh"]
        @test ax1.limits[][2] == (0, 31) # limits set on ax1, not on the current axis
    end

    @testset "vecplot" begin
        P = IF(N = 3, name = "P")
        P.I .= 400pA
        monitor!(P, [:v, :w])
        sim!([P]; duration = 100ms)
        fig, ax, plt = vecplot(P, [:v])
        @test plt.label[] == "v"
        @test vecplot([P, P], :v) isa Figure
        @test vecplot(P, Any[:v, :w]) isa Figure            # one panel per symbol
        @test vecplot(P, [:v, :w]) isa Makie.FigureAxisPlot  # all symbols in one axis
        f = Figure(); ax = Axis(f[1, 1])
        @test vecplot!(ax, P, :v; factor = 2.0f0, neurons = 1:2) !== nothing
        @test vecplot!(ax, P, :v; factor = :v, neurons = 1:2) !== nothing
    end
end
