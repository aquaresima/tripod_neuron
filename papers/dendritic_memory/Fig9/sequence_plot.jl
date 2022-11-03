"This file reproduces Fig 10c and 10d"

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include("stimuli.jl")
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

using StatsBase, JLD2, LaTeXStrings
##
# file = h5open(datadir("dendritic_memory","sequence_scores.h5")) do file
file = datadir("dendritic_memory","sequence_data.jld2")
JLD2.load(file)
scores = JLD2.load(file)["data"]
inputs = JLD2.load(file)["inputs"]
@unpack βs, intervals = inputs


plots =[]
freqs = reverse(1000 ./intervals)
for (scores_values, label) in zip(values(scores), keys(scores))
    # @show model, label
        p = heatmap(freqs, βs , reverse(scores_values, dims=1)', label=label, clims=(0.5,1), xscale=:log, c=:amp)
        length(plots) == 2 &&(plot!(ylabel="Inhibition ratio (μ)", xlabel="Switch frequency (Hz) "))
        plot!(title=label)
        push!(plots,p)
end
p =plot(plots..., cbar=false, size=(600,600))

savefig(p, plotsdir("dendritic_memory","Fig9G.pdf"))