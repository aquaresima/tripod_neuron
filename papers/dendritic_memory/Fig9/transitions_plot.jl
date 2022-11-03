"""
This file reproduces Fig 9C and 9D
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include("stimuli.jl")
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using JLD2

file =datadir("dendritic_memory","transitions_corr_data_#1.jld2")

fid = JLD2.load(file)
scores_oscillation = fid["scores_oscillation"]
scores_signal = fid["scores_exp_decay"]
intervals = fid["intervals"]
shifts =fid["shifts"]

oscillation = maximum(abs.(scores_oscillation),dims=3)[:,:,1]
signal = maximum(abs.(scores_signal),dims=3)[:,:,1]

frequencies = 1000 ./intervals

c = collect(palette(:roma,4)) |> x-> begin x[end] = RGBA(0,0,0,0); return x' end



p = plot(plot(frequencies, oscillation', legend=false,  marker=:circle,msc=:transparent, ms=5, xscale=:log, ylims=(0,1),c=c, xlabel="                            Switch frequency (ms)",),
    plot(frequencies,signal', marker=:circle, ms=5, msc=:transparent, xscale=:log, ylims=(0,1), c=c),
    plot(shifts,scores_oscillation[:,15,:]', ylims=(-0.9,1), c=c, ylabel="                       Signal/spikes correlation", xlabel="                           Response delay (ms)"),
    plot(shifts,scores_signal[:,15,:]', ylims=(-0.9,1.), c=c),
    legend=false, lw=2, margin=5mm
    )

savefig(p, plotsdir("dendritic_memory","Fig9CD.pdf"))
p
