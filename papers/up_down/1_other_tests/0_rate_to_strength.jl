using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsPlots, JLD2
include(projectdir("scripts", "up_down", "default_plots.jl"))

function make_spikes_bursty(simtime, d1, d2, f; nmda, γ = 0.1)
    # @assert γ>=0 && γ<=1
    cond = TN.get_balance_conditions(d1, d2, f, nmda = nmda)
    kie = TN.get_balance_conditions(d1, d2, f, nmda = nmda).ie_ratio
    istdp_syn =
        TN.get_balance_conditions(d1, d2, f, nmda = nmda, istdp_syn = true).synapses[end]
    # @show kie, istdp_syn
    istdp_syn = 1 - (1 - istdp_syn) * (1 - γ)
    kie = 1 .- (1 .- kie) * γ
    new_cond = (cond.rate, ie_ratio = kie)
    inputs = TN.make_spikes(simtime; new_cond...)
    synapses = [1, 1, 1, 1, istdp_syn, istdp_syn]
    # @show kie, synapses
    # @show kie[end]*synapses[end]
    # @show kie[end] * synapses[end]
    return inputs, synapses
end



γs = 0:0.05:1
simtime = 11_000
for nmda in [true, false]
    mean_v = zeros(3, length(γs), length(TN.ds), length(TN.νs))
    std_v = zeros(3, length(γs), length(TN.ds), length(TN.νs))
    rate = zeros(2, length(γs), length(TN.ds), length(TN.νs))

    for f in eachindex(TN.νs)
        @info "Rate:", TN.νs[f]
        Threads.@threads for d in eachindex(TN.ds)
            for g in eachindex(γs)
                γ = γs[g]
                # rate balance	
                model =
                    (ds = (TN.ds[d], TN.ds[d]), nmda = nmda, soma_only = sum(TN.ds[d]) == 0)
                inputs, synapses = make_spikes_bursty(simtime, d, d, f, nmda = true, γ = γ)
                v_spikes = TN.run_tripod(
                    inputs,
                    simtime;
                    model...,
                    do_spikes = true,
                    synapses = synapses,
                )
                v_nospikes = TN.run_tripod(
                    inputs,
                    simtime;
                    model...,
                    do_spikes = false,
                    synapses = synapses,
                )
                mean_v[:, g, d, f] = mean(v_nospikes[:, 100_000:end], dims = 2)[:, 1]
                std_v[:, g, d, f] = std(v_nospikes[:, 100_000:end], dims = 2)[:, 1]
                rate[1, g, d, f] = TN.get_spike_rate(v_spikes)
                rate[2, g, d, f] = TN.get_cv(v_spikes)
            end
        end
    end
    syn = nmda ? "NMDA" : "AMPA"
    data = @strdict mean_v std_v rate
    filep = datadir("up_down", "activity", "rate_to_strength_$syn.jld2")
    safesave(filep, data)
end


# inputs, synapses = make_spikes_bursty(1000,300,300,34,nmda=true, γ=0)
