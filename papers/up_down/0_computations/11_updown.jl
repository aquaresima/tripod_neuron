using DrWatson
@quickactivate "Tripod"

using TripodNeuron
using JLD2
using ProgressBars, Logging

include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
include(projectdir("scripts", "plots", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_analysis.jl"))

##
@unpack νs = TN
simtime = 100_000
βs = 0:50:500
dend_model = TN.models[2]
interval = 150_000:200_000
data = Array{Any,3}(undef, 3, length(βs), length(νs))
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
for n in ProgressBar(eachindex(νs))
    Threads.@threads for m in eachindex(βs)

        model = dend_model
        cond = TN.get_balance_conditions(model.ds..., n; nmda = true)
        inputs = TN.make_spikes(simtime, βs[m]; cond...)
        voltage =
            TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
        data[1, m, n] = updown_analysis(voltage[1, :])

        model = dend_model
        cond = TN.get_balance_conditions(model.ds..., n; nmda = false)
        inputs = TN.make_spikes(simtime, βs[m]; cond...)
        voltage =
            TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
        data[2, m, n] = updown_analysis(voltage[1, :])

        model = TN.models[4]
        cond = TN.get_balance_conditions(model.ds..., n; nmda = false)
        inputs = TN.make_spikes(simtime, βs[m]; cond...)
        voltage =
            TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
        data[3, m, n] = updown_analysis(voltage[1, :])
    end
end

data = @strdict data βs νs model = dend_model
save(datadir("up_down", "bistability", "characterize_updown_$(dend_model.ds).jld2"), data)
