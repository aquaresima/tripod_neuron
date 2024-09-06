using DrWatson
using Plots, Revise, Random
@quickactivate "Tripod"
using TripodNeuron
using Statistics

##
β = 200
simtime = 25_000
Random.seed!(121)
do_spikes = false
model = (350, 350, 13)
cond = TN.get_balance_conditions(model..., nmda = true)
inputs = TN.make_spikes(simtime, β; cond...)
voltage = TN.run_tripod(
    cond.ds,
    inputs,
    simtime,
    do_spikes = do_spikes,
    soma_only = cond.soma_only,
    syn_model = cond.syn_model,
)
output = round(TN.get_spike_rate(voltage[1, :]), digits = 2)
μ = mean(voltage[:, 5000:end], dims = 2)[:, 1]
plot(
    [
        plot(voltage[n, :]) |> p -> hline!(p, μ[n:n], ylims = (-70, -10), legend = false)
        for n = 1:3
    ]...,
    layout = (3, 1),
)
title = "νₒᵤₜ: $output Hz, μ: $(round(μ[1], digits=2)) mV "
p1 = plot!(subplot = 1, title = title, ylims = (-70, -20))
p1 = plot!(xticks = :none)
p1 = plot!(
    subplot = 3,
    xticks = (
        range(1, simtime * 10, length = 5),
        round.(range(1, simtime / 1000, length = 5), digits = 2),
    ),
    xlabel = "Time (s)",
)
# plot!(ylims=(50_000, 70_000))
plot!(size = (600, 800))
histogram(voltage[1, :], bins = 100, title = "Soma voltage")

# cor(voltage[3,:] , voltage[2,:])
