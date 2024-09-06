using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
using JLD2
using ProgressBars, Logging


##

## From critical window analysis

function get_balance_window(simtime, input_rates, ratios; model, nmda, β = 200, samples = 5)
    c_window = zeros(samples, 4, length(input_rates), length(ratios))
    μ = zeros(samples, 3, length(input_rates), length(ratios))
    output_rates = zeros(samples, length(input_rates), length(ratios))
    AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)

    for n in ProgressBar(eachindex(input_rates))
        Threads.@threads for m in eachindex(ratios)
            cond = TN.get_balance_conditions(model.ds..., n; nmda = nmda, shift = ratios[m])
            for z = 1:samples
                inputs = TN.make_spikes(simtime, β; cond...)
                voltage = TN.run_tripod(
                    inputs,
                    simtime;
                    model...,
                    AdEx = AdEx,
                    syn_model = cond.syn_model,
                )
                c0 = critical_window(voltage[1, :], ratio = 0.0)
                c1 = critical_window(voltage[1, :], ratio = 0.10)
                c3 = critical_window(voltage[1, :], ratio = 0.30)
                c5 = critical_window(voltage[1, :], ratio = 0.50)
                #             all_w = all_windows(voltage[1,:], max_b = max_b)
                output_rates[z, n, m] = round(TN.get_spike_rate(voltage[1, :]), digits = 2)
                μ[z, :, n, m] = mean(voltage[:, 5000:end], dims = 2)[:, 1]
                c_window[z, 1, n, m] = c0
                c_window[z, 2, n, m] = c1
                c_window[z, 3, n, m] = c3
                c_window[z, 4, n, m] = c5
            end
        end
    end
    return (
        critical = mean(c_window, dims = 1)[1, :, :, :],
        output_rates = mean(output_rates, dims = 1)[1, :, :, :],
        μ = mean(μ, dims = 1)[1, :, :, :],
        β = β,
    )
end


ratios = collect(0.04:0.04:2)
input_rates = TN.νs
simtime = 50_000
β = 200

@info length(TN.νs) TN.μmem

#nmda =true
for nmda in [true, false]
    syn = nmda ? "NMDA" : "AMPA"
    data = Dict()
    for (model, label) in zip(TN.models, TN.labels)
        @info "Model: $(label) with $syn"
        data_model = get_balance_window(
            simtime,
            input_rates,
            ratios;
            model = model,
            nmda = nmda,
            β = β,
        )
        push!(data, model => data_model)
    end
    data = @strdict data ratios input_rates
    file = datadir("up_down", "bistability", "balance_$syn.jld2")
    safesave(file, data)
end
