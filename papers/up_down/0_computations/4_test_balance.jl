using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsPlots, JLD2

simtime = 50_000
νs = TN.νs
for nmda in [true, false]
    mean_v = zeros(3, length(TN.ds), length(TN.ds), length(νs))
    std_v = zeros(3, length(TN.ds), length(TN.ds), length(νs))
    rate = zeros(2, length(TN.ds), length(TN.ds), length(νs))
    istdp_mean_v = zeros(3, length(TN.ds), length(TN.ds), length(νs))
    istdp_std_v = zeros(3, length(TN.ds), length(TN.ds), length(νs))
    istdp_rate = zeros(2, length(TN.ds), length(TN.ds), length(νs))

    for f in eachindex(νs)
        @info "Rate:", νs[f]
        Threads.@threads for d1 in eachindex(TN.ds)
            for d2 in eachindex(TN.ds)
                if d1 == 1 || d2 == 1
                    model = (
                        ds = (TN.ds[1], TN.ds[1]),
                        nmda = nmda,
                        soma_only = sum(TN.ds[1]) == 0,
                    )
                end
                # rate balance
                model = (
                    ds = (TN.ds[d1], TN.ds[d2]),
                    nmda = nmda,
                    soma_only = sum(TN.ds[d1]) == 0,
                )
                cond = TN.get_balance_conditions(d1, d2, f, nmda = nmda)
                inputs = TN.make_spikes(simtime; cond...)
                v_spikes =
                    TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
                v_nospikes =
                    TN.run_tripod(inputs, simtime; model..., do_spikes = false, cond...)
                mean_v[:, d1, d2, f] = mean(v_nospikes[:, 200_000:end], dims = 2)[:, 1]
                std_v[:, d1, d2, f] = std(v_nospikes[:, 200_000:end], dims = 2)[:, 1]
                rate[1, d1, d2, f] = TN.get_spike_rate(v_spikes)
                rate[2, d1, d2, f] = TN.get_cv(v_spikes)
                ## istdp balance
                cond = TN.get_balance_conditions(d1, d2, f, nmda = nmda, istdp_syn = true)
                inputs = TN.make_spikes(simtime; cond...)
                v_spikes =
                    TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
                v_nospikes =
                    TN.run_tripod(inputs, simtime; model..., do_spikes = false, cond...)
                istdp_mean_v[:, d1, d2, f] =
                    mean(v_nospikes[:, 200_000:end], dims = 2)[:, 1]
                istdp_std_v[:, d1, d2, f] = std(v_nospikes[:, 200_000:end], dims = 2)[:, 1]
                istdp_rate[1, d1, d2, f] = TN.get_spike_rate(v_spikes)
                istdp_rate[2, d1, d2, f] = TN.get_cv(v_spikes)
                @warn "Mean: ", d1, d2, f, mean_v[1, d1, d2, f]
            end
        end
    end
    syn = nmda ? "NMDA" : "AMPA"
    data = @strdict mean_v std_v rate istdp_mean_v istdp_std_v istdp_rate
    filep = datadir("up_down", "activity", "full_activity_$syn.jld2")
    safesave(filep, data)
end
