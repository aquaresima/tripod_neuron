using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsPlots, JLD2, ProgressBars, Printf
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))

simtime = 50_000
νs = TN.νs

for β in [200, 400]
    for nmda in [true, false]
        rate = zeros(2, length(TN.ds), length(TN.ds), length(νs))
        cors = zeros(3, length(TN.ds), length(TN.ds), length(νs))
        c_window = zeros(2, length(TN.ds), length(TN.ds), length(νs))
        AdEx = TN.AdExParams(idle = 0.1, up = 0.1)

        iter = ProgressBar(eachindex(νs))
        for f in iter
            @info "Rate:", νs[f]
            Threads.@threads for d2 in eachindex(TN.ds)
                for d1 in eachindex(TN.ds)
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
                    inputs = TN.make_spikes(simtime, β; cond...)

                    v_spikes = TN.run_tripod(
                        inputs,
                        simtime;
                        model...,
                        do_spikes = true,
                        cond...,
                        AdEx = AdEx,
                    )

                    v_spikes = v_spikes[:, 50_000:end]
                    rate[1, d1, d2, f] = TN.get_spike_rate(v_spikes)
                    rate[2, d1, d2, f] = TN.get_cv(v_spikes)

                    c1 = cor(v_spikes[2, :], v_spikes[3, :])
                    c2 = cor(v_spikes[1, :], v_spikes[2, :])
                    c3 = cor(v_spikes[1, :], v_spikes[3, :])
                    cors[:, d1, d2, f] = [c1, c2, c3]

                    c1 = critical_window(v_spikes[1, :], ratio = 0.10)
                    c2 = critical_window(v_spikes[1, :], ratio = 0.30)
                    c_window[1, d1, d2, f] = c1
                    c_window[2, d1, d2, f] = c2
                    set_description(
                        iter,
                        string(
                            @sprintf(
                                "Dendrites: %d %d %d Mean: %.2f ",
                                d1,
                                d2,
                                f,
                                rate[1, d1, d2, f]
                            )
                        ),
                    )
                end
            end
        end
        syn = nmda ? "NMDA" : "AMPA"
        data = @strdict rate β nmda cors c_window
        filep = datadir("up_down", "asymmetries", "rate_β=$(β)$(syn).jld2")
        safesave(filep, data)
    end
end
