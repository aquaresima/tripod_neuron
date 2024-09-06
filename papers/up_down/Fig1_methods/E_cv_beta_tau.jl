using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsBase, StatsPlots, JLD2, Plots, RollingFunctions
include(projectdir("scripts", "plots", "default_plots.jl"));
using MultivariateStats
using Plots

##

data = collect(values(load(datadir("up_down", "inputs", "cv_inputs.bson"))))[1]
@unpack cvs, τs, βs = data
# plot(τs[2:end], cvs[2:end,:], c=my_palette', label="")
##
# plot(βs ./1000, cvs[:,:]', c=my_palette', label="")
q2 = contour(
    βs[2:end],
    τs,
    cvs[:, 2:end, 1],
    title = "CV of input spike trains",
    xlabel = "Fluctuations (β)",
    ylabel = "τ (ms)",
    legend = :topleft,
    yscale = :log10,
    xscale = :log10,
    fill = true,
    c = :amp,
    rightmargin = 10mm,
)
plot!(ylabel = "CV", xlabel = "β/1000", legend = :topleft, title = "CV of input spikes")
# heatmap(τs, my_palette, title="τ (ms)", label="τ (ms)", legend=:topleft)

plot(TN.make_rates(10000, 10; τ = 50))
plot!(TN.make_rates(10000, 10; τ = 50))
plot!(TN.make_rates(10000, 10; τ = 50))
mean([std(TN.make_rates(10000, 10; τ = 50)) for _ = 1:10])
variance = [mean([std(TN.make_rates(10000, β; τ = 500)) for _ = 1:10]) for β = 0:50:1000]
plot(0:50:1000, variance)

plot(1:1000, yscale = :log, xscale = :log)
