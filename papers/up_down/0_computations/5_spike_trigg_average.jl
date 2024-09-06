using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, Logging, JLD2
##


spiketime = 5_000
samples = 500
function get_STA(model)
    if model.nmda
        local syn_model = TN.human_synapses
    else
        local syn_model = TN.ampa_equivalent
    end
    simtime = 1_000
    voltage = TN.run_tripod((0, 0), zeros(6), simtime, soma_only = true)
    volt = zeros(size(voltage)..., length(TN.νs), samples)
    Threads.@threads for n in eachindex(TN.νs)
        for s = 1:samples
            cond = TN.get_balance_conditions(model.ds..., n; nmda = model.nmda)
            local inputs = TN.make_spikes(simtime; cond...)
            local voltage = TN.run_tripod(
                inputs,
                simtime;
                model...,
                syn_model = syn_model,
                do_spikes = false,
            )
            local inputs[2, spiketime] += 2.0
            local voltage_1 = TN.run_tripod(
                inputs,
                simtime;
                model...,
                syn_model = syn_model,
                do_spikes = false,
            )
            volt[:, :, n, s] .= voltage_1 - voltage
        end
    end
    return mean(volt, dims = 4)[:, :, :, 1]
end

function get_peak_and_width(model)
    voltage = get_STA(model)
    peak = maximum(voltage[1, :, :], dims = 1)
    half = [_half_width(v[:, 1]) for v in eachslice(voltage[1, :, :], dims = 2)]
    return (mem = voltage, peak = peak, half = half)
end

function _half_width(voltage)
    f = argmax(voltage)
    m = maximum(voltage)
    l = findfirst(x -> x < m / 2, voltage[f:end])
    l = isnothing(l) ? length(voltage) : l
    return (l + f) * 0.1
end


##
@info "Run Spike triggered average for short-prox models"
soma_NMDA = (ds = (0, 0), nmda = true, soma_only = true) |> get_peak_and_width
soma_AMPA = (ds = (0, 0), nmda = false, soma_only = true) |> get_peak_and_width
dist_NMDA = (ds = (400, 400), nmda = true) |> get_peak_and_width
prox_NMDA = (ds = (150, 150), nmda = true) |> get_peak_and_width
dist_AMPA = (ds = (400, 400), nmda = false) |> get_peak_and_width
prox_AMPA = (ds = (150, 150), nmda = false) |> get_peak_and_width
data = @strdict spiketime dist_NMDA prox_NMDA dist_AMPA prox_AMPA soma_AMPA soma_NMDA
filep = datadir("up_down", "activity", "added_spike.jld2")
safesave(filep, data)
