##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
using EllipsisNotation
using Plots
using Random


function alpha_syn(t; τ, A)
    return maximum([A * (t - 20) * exp(1 - (t - 20) / τ), 0])
end

τ = 0.2
A = 7.1 / τ

plot(
    alpha_syn.(0:0.1:60; τ = τ, A = A),
    label = " α Kuhn synapse",
    xlabel = "Time (ms)",
    ylabel = "Syn. conductance (nS)",
    legend = :bottomright,
)

# @eval TN begin
#     EsynSoma∏

Random.seed!(1234)
simtime = 60
model = TN.models[4]
cond = TN.get_balance_conditions(model.ds..., 20; syn_model = :kuhn)
# inputs = 
# keys(cond)
inputs = zeros(Int, 6, 1000)
inputs[2, 200] = 1
voltage, syn =
    TN.run_tripod(model.ds, inputs, simtime; do_spikes = true, cond..., rec_syn = true)
sum(syn)
plot!(syn[2, :, 1], label = " double exp \n Kuhn synapse")
plot!(legend = :topright, size = (600, 400))
