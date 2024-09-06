using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using JLD, Statistics

##


νs, ds = TN.νs, TN.ds
AdEx = TN.AdExParams(u_r = -55)

simtime = 90_000
rec_ws = zeros(3, length(νs), length(ds))
rec_rate = zeros(length(νs), length(ds))
Threads.@threads for n in eachindex(νs)
    @show νs[n]
    for m in eachindex(ds)
        ν = νs[n]
        simtime = round(Int, 50_000 + (50_000 / ν))
        d = ds[m]
        inputs = [0.0, ν, ν, 0.0, ν, ν]
        v, w = TN.run_tripod_vogels(
            (d, d),
            inputs,
            simtime,
            AdEx = AdEx,
            soma_only = (sum(d) == 0),
        )
        rec_rate[n, m] = TN.get_spike_rate(v[1, end-10_000:end])
        rec_ws[:, n, m] .= mean(w[:, end-10_000:end], dims = 2)[:, 1]
    end
end

file = datadir("up_down", "balance", "inhibitory_balance_istdp_NMDA_ur$(AdEx.u_r).jld2")
safesave(file, @strdict data = rec_ws rate = rec_rate)

rec_ws = zeros(3, length(νs), length(ds))
rec_rate = zeros(length(νs), length(ds))
Threads.@threads for n in eachindex(νs)
    @show νs[n]
    for m in eachindex(ds)
        ν = νs[n]
        simtime = round(Int, 50_000 + (50_000 / ν))

        d = ds[m]
        inputs = [0.0, ν, ν, 0.0, ν, ν]
        v, w = TN.run_tripod_vogels(
            (d, d),
            inputs,
            simtime,
            syn_model = TN.ampa_equivalent,
            AdEx = AdEx,
            soma_only = (sum(d) == 0),
        )
        rec_rate[n, m] = TN.get_spike_rate(v[1, end-10_000:end])
        rec_ws[:, n, m] .= mean(w[:, end-10_000:end], dims = 2)[:, 1]
    end
end

file = datadir("up_down", "balance", "inhibitory_balance_istdp_AMPA_ur$(AdEx.u_r).jld2")
safesave(file, @strdict data = rec_ws rate = rec_rate)
