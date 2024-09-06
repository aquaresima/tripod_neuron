using DrWatson
@quickactivate "Tripod"

using TripodNeuron
using JLD2
using ProgressBars, Logging
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
include("inhibitory_balance.jl")


##
simtime = 100_000
βs = 0:50:500
model = TN.models[3]
interval = 150_000:200_000
NMDAs = Array{Any,2}(undef, length(βs), length(νs))
Threads.@threads for n in eachindex(νs)
    for m in eachindex(βs)
        cond = TN.get_balance_conditions(model.ds..., n; nmda = true)
        inputs = TN.make_spikes(simtime, βs[m]; cond...)
        AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
        voltage =
            TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
        NMDA = updown_analysis(voltage[1, :])
        NMDAs[m, n] = NMDA
    end
end
