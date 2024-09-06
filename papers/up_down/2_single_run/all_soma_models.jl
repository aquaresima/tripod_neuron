using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsPlots, JLD2, Plots, RollingFunctions
include(projectdir("scripts", "up_down", "default_plots.jl"))
using Random
using Plots

##
Random.seed!(1234)
simtime = 1_0000
model = TN.models[4]
_rate = 5

plot()
Random.seed!(1234)
for syn in [TN.nmda_soma, TN.ampa_kuhn, TN.human_synapses, TN.ampa_equivalent]
    cond = TN.get_balance_conditions(model.ds..., _rate; nmda = false, syn_model = syn)
    inputs = TN.make_spikes(simtime; cond...)
    Kie = 0.30
    ν = TN.νs[_rate]
    inputs = [0.0, ν, ν, 0.0, Kie * ν, Kie * ν]
    voltage = TN.run_tripod(
        model.ds,
        inputs,
        simtime;
        do_spikes = true,
        synapses = TN.ampa_equivalent,
        cond...,
    )
    print(sum(inputs, dims = 2))
    plot!(voltage[1, :], title = "Soma models", label = string(syn.name))
end
plot!()
##
