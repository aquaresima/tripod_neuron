"""
Compute inhibitory balance for the tripod model.
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars


νs = collect(exp.(range(log(0.5), log(80), length = 44)))
@info "min ν: $(round(minimum(νs),digits=2)), max ν: $(round(maximum(νs),digits=2)) length ν: $(length(νs))"

function compute_voltage_test(; model, μ_mem::Real, ν::Real, Kie::Float64, v::Bool = false)
    simtime = 5000
    inputs = [0.0, ν, ν, 0.0, Kie * ν, Kie * ν]
    voltage = TN.run_tripod(inputs, simtime; model...)
    μ_computed = mean(voltage[1, 10000:end])
    (v) && (return voltage)
    target = abs(μ_mem - μ_computed)
    return target
end

function compute_voltage(; model, μ_mem::Real, ν::Real, Kie::Float64, v::Bool = false)
    simtime = 10000
    inputs = [0.0, ν, ν, 0.0, Kie * ν, Kie * ν]
    voltage = TN.run_tripod(inputs, simtime; model...)
    μ_computed = mean(voltage[1, 10000:end])
    (v) && (return voltage)
    target = abs(μ_mem - μ_computed)
    return target
end

function compute_balance_condition(
    μmem,
    nmda = true;
    ds = TN.ds,
    νs = νs,
    syn_model::TN.TripodSynapses,
)
    kies = collect(0.00:0.01:4.00)
    if isnothing(syn_model)
        syn_model = nmda ? TN.human_synapses : TN.ampa_equivalent
    else
        syn_model = syn_model
    end
    data = zeros(length(νs), length(ds), length(kies))
    minμ = 10
    for m in ProgressBar(eachindex(ds))
        Threads.@threads for i in eachindex(νs)
            model = (
                ds = (ds[m], ds[m]),
                soma_only = sum(ds[m]) == 0 ? true : false,
                syn_model = syn_model,
                AdEx = TN.AdEx,
                do_spikes = false,
            )

            # μmem = model.soma_only ? -55. : -50.
            # μmem = 
            allμ = minμ
            for k in eachindex(kies)
                μ = compute_voltage_test(
                    μ_mem = μmem,
                    model = model,
                    ν = νs[i],
                    Kie = kies[k],
                )
                if μ < minμ
                    for _ = 1:5
                        μ = min(
                            compute_voltage(
                                μ_mem = μmem,
                                model = model,
                                ν = νs[i],
                                Kie = kies[k],
                            ),
                            μ,
                        )
                        (μ < 0.5) && (break)
                    end
                    data[i, m, k] = μ
                    allμ = min(μ, allμ)
                else
                    data[i, m, k] = minμ
                end
            end
            @show round(νs[i], digits = 1), ds[m], allμ
        end
    end
    return data, νs, ds, kies
end
