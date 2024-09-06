using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
using JLD2
using ProgressBars, Logging


##
function get_critical_window(simtime, input_rates, βs; model, nmda, samples = 5)
    c_window = zeros(samples, 4, length(input_rates), length(βs))
    μ = zeros(samples, 3, length(input_rates), length(βs))
    output_rates = zeros(samples, length(input_rates), length(βs))
    AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)

    for n in ProgressBar(eachindex(input_rates))
        cond = TN.get_balance_conditions(model.ds..., n, nmda = nmda)
        Threads.@threads for m in eachindex(βs)
            for z = 1:samples
                β = βs[m]
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
    )
end


βs = 0:20:1000
input_rates = TN.νs
simtime = 50_000
for nmda in [true, false]
    data = Dict()
    syn = nmda ? "NMDA" : "AMPA"
    for (model, label) in zip(TN.models, TN.labels)
        @info "Model: $(label) with $syn"
        data_model =
            get_critical_window(simtime, input_rates, βs; model = model, nmda = nmda)
        push!(data, label => data_model)
    end
    data = @strdict data βs = βs
    file = datadir("up_down", "bistability", "critical_window_$syn.jld2")
    safesave(file, data)
end


#data_ampa = get_bimodal_window(βs, input_rates, simtime; nmda=false)
#safesave(datadir("up_down","robustness","critical_window_ampa_largeβ.bson"),@dict data_ampa)
