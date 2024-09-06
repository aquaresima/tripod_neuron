using DrWatson
@quickactivate "Tripod"
using TripodNeuron

using Plots, Revise, Random
using Statistics
include(projectdir("scripts", "up_down", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_plots.jl"));
include(projectdir("scripts", "up_down", "updown_analysis.jl"))
##



##
# NMDA = updown_analysis(voltage[1, :])

simtime = 50_000
βs = 100:20:500
@unpack νs = TN


AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
for nmda in [true, false]
    syn = nmda ? "NMDA" : "AMPA"
    data = Dict()
    for (model, label) in zip(TN.models, TN.labels)
        @info "$label $syn $model"
        # Threads.@threads 
        local states = zeros(length(νs), length(βs), 2, 4)
        for f in eachindex(TN.νs)
            for b in eachindex(βs)
                β = βs[b]
                cond = TN.get_balance_conditions(model.ds..., f; nmda = nmda)
                inputs = TN.make_spikes(simtime, β; cond...)
                voltage = TN.run_tripod(
                    inputs,
                    simtime;
                    model...,
                    cond...,
                    AdEx = AdEx,
                    do_spikes = true,
                )
                state = updown_analysis(voltage[1, :])
                local states[f, b, 1, :] .= mean(state.up, dims = 2)[:, 1]
                local states[f, b, 2, :] .= mean(state.down, dims = 2)[:, 1]
            end
        end
        push!(data, label => states)
    end
    data = @strdict data βs
    file = datadir("up_down", "robustness", "updown_$syn.jld2")
    safesave(file, data)
end
# #
