using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
using JLD2
using ProgressBars, Logging


##
function get_critical_window(simtime, τs, βs, input_rate; model, nmda, samples = 5)
    c_window = zeros(samples, 2, length(τs), length(βs))
    μ = zeros(samples, 3, length(τs), length(βs))
    output_rates = zeros(samples, length(τs), length(βs))
    AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
    autocor_times = 1:10:1000
    correlations = zeros(samples, length(τs), length(βs), 3, length(autocor_times))
    crosscorr = zeros(samples, length(τs), length(βs), 3)

    for n in ProgressBar(eachindex(τs))
        cond = TN.get_balance_conditions(model.ds..., input_rate, nmda = nmda)
        Threads.@threads for m in eachindex(βs)
            for z = 1:samples
                β = βs[m]
                inputs = TN.make_spikes(simtime, β; τ = τs[n], cond...)
                voltage = TN.run_tripod(
                    inputs,
                    simtime;
                    model...,
                    AdEx = AdEx,
                    syn_model = cond.syn_model,
                )
                c1 = critical_window(voltage[1, :], ratio = 0.10)
                c3 = critical_window(voltage[1, :], ratio = 0.30)
                c_window[z, 1, n, m] = c1
                c_window[z, 2, n, m] = c3

                output_rates[z, n, m] = round(TN.get_spike_rate(voltage[1, :]), digits = 2)

                μ[z, :, n, m] = mean(voltage[:, 5000:end], dims = 2)[:, 1]

                correlations[z, n, m, 1, :] = autocor(voltage[1, 5000:end], autocor_times)
                correlations[z, n, m, 2, :] = autocor(voltage[2, 5000:end], autocor_times)
                correlations[z, n, m, 3, :] = autocor(voltage[3, 5000:end], autocor_times)

                c1 = cor(voltage[2, :], voltage[3, :])
                c2 = cor(voltage[1, :], voltage[2, :])
                c3 = cor(voltage[1, :], voltage[3, :])
                crosscorr[z, n, m, :] = [c1, c2, c3]
            end
        end
    end
    return (
        critical = mean(c_window, dims = 1)[1, :, :, :],
        output_rates = mean(output_rates, dims = 1)[1, :, :, :],
        μ = mean(μ, dims = 1)[1, :, :, :],
        autocor = mean(autocor_times, dims = 1)[1, :, :, :, :],
        crosscorr = mean(crosscorr, dims = 1)[1, :, :, :, :],
    )
end


begin
    βs = 0:20:1000
    τs = collect(exp.(range(log(0.1), log(1000), length = 50)))
    @unpack νs = TN
    input_rate = 6#; νs[input_rate]
    simtime = 50_000
    for input_rate in [12, 18, 6, 24]
        for nmda in [true]
            data = Dict()
            syn = nmda ? "NMDA" : "AMPA"
            for (model, label) in zip(TN.models, TN.labels)
                @info "Model: $(label) with $syn"
                data_model = get_critical_window(
                    simtime,
                    τs,
                    βs,
                    input_rate;
                    model = model,
                    nmda = nmda,
                )
                push!(data, label => data_model)
            end
            data = @strdict data βs = βs τs = τs ν = νs[6]
            file = datadir(
                "up_down",
                "bistability",
                "critical_beta_tau_$(syn)_rate=$input_rate.jld2",
            )
            safesave(file, data)
        end
    end
end


#data_ampa = get_bimodal_window(βs, τs, simtime; nmda=false)
#safesave(datadir("up_down","robustness","critical_window_ampa_largeβ.bson"),@dict data_ampa)
