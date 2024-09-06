"""
Measure the PSP_peak by comparing the extreme value reached by the `voltage` and the `rest` (70.6 mV) potential. By default consider the maximum depolarized value, hence the EPSP.

# Arguments:
spiketime : the first time-index to consider for the EPSP peak.
rest : the resting potential of the neuron.
IPSP : measure neuron hyperpolarization instead of depolarization (IPSP).
"""

function get_PSP_peak(
    voltage::Array{Float32,2};
    spiketime = -1,
    rest = -70.6,
    IPSP = false,
    compartment = 1,
)
    spiketime = spiketime < 0 ? EXCSPIKETIME : spiketime
    if IPSP == false
        return maximum(voltage[compartment, spiketime:end]) - rest
    else
        return minimum(voltage[compartment, spiketime:end]) - rest
    end
end

function _PoissonInput(Hz_rate::Real, interval::Int64, dt::Float32)
    λ = 1000 / Hz_rate
    spikes = zeros(Int, round(Int, interval / dt))
    t = 1
    while t < interval / dt
        Δ = rand(Exponential(λ / dt))
        t += Δ
        if t < interval / dt
            spikes[round(Int, t)] += 1
        end
    end
    return spikes
end

function PoissonInput(Hz_rate::Real, interval::Int64, dt::Float32; neurons::Int64 = 1)
    spikes = zeros(Int, neurons, round(Int, interval / dt))
    for n = 1:neurons
        spikes[n, :] .= _PoissonInput(Hz_rate::Real, interval::Int64, dt::Float32)
    end
    return spikes
end


function null_input(; kwargs...) end
function null_input(a, b; kwargs...)
    return 0
end


function make_rates(simtime::Int, β = 0, dt = 0.1; rate::Real = 1000.0, τ = 50, kwargs...)
    """
    Produce inpute rates, dimension is kHz⁻¹
    Resolution is ms.
    """
    @debug "Generating input rates with τ = $(round(τ, digits=1)) ms and β = $β, rate = $(round(rate, digits=1)) Hz"
    cc = exp(-dt / τ)
    noise = 0.0
    inputs = zeros(round(Int, simtime))
    for t = 1:round(Int, simtime)
        re = rand() - 0.5
        noise = (noise - re) * cc + re
        inputs[t] = 1 + maximum([0, noise]) * β
    end
    return inputs / (sum(inputs) / simtime) * rate
end

function make_spikes(
    simtime,
    β = 0,
    ;
    rate = 1000.0,
    ie_ratio = ones(3),
    rates = nothing,
    seed = nothing,
    kwargs...,
)
    r_exc = zeros(3, simtime)
    r_inh = zeros(3, simtime)
    if isnothing(rates)
        for i = 1:3
            r_exc[i, :] = make_rates(simtime, β, rate = rate; kwargs...)
            r_inh[i, :] = make_rates(simtime, 0, rate = rate; kwargs...)
        end
    else
        r_exc, r_inh = rates
    end
    if !isnothing(seed)
        Random.seed!(seed)
    end
    exc_spikes = zeros(Float32, 3, round(Int, simtime / dt))
    inh_spikes = zeros(Float32, 3, round(Int, simtime / dt))
    for x = 1:simtime
        for y = 1:round(Int, 1 / dt)
            z = (x - 1) * 10 + y
            for i = 2:3
                exc_spikes[i, z] = rand(Poisson(r_exc[i, x] * dt / 1000))
                inh_spikes[i, z] = rand(Poisson(ie_ratio[i] * r_inh[i, x] * dt / 1000))
            end
        end
    end
    return vcat(exc_spikes, inh_spikes)
end

function make_balance(
    simtime;
    rate = 3000,
    ie_ratio = [0.0, 1.0, 1],
    model = nothing,
    kwargs...,
)
    isnothing(model) || (ie_ratio = get_ie_ratio(model))
    return make_spikes(simtime, β = 0, ; rate = rate, ie_ratio = ie_ratio, kwargs...)
end

# Balanced condition:
## With this rate and ie_ratio for the dendrites the soma rests around -50 mV and have a low firing rate.
function get_ie_ratio(model)
    ie_long = 1.0
    ie_short = 1.6
    ie_soma = 0.3
    ie_ratio = [0.0]
    for d in model.ds
        if d == 400
            push!(ie_ratio, ie_long)
        elseif d == 150
            push!(ie_ratio, ie_short)
        elseif d == 0
            push!(ie_ratio, ie_soma)
        else
            throw("The set dendritic length has no balance ratio")
        end
    end
    return ie_ratio
end

#
get_spikes(voltage::Array{T,1}) where {T<:Real} =
    [v == postspike.AP_membrane for v in voltage[:]]
get_spike_times(voltage::Vector{T}) where {T<:Real} = findall(get_spikes(voltage))
get_spike_times(voltage::Matrix{T}) where {T<:Real} = get_spike_times(voltage[1, :])
get_spike_rate(voltage::Vector{T}) where {T<:Real} =
    length(get_spike_times(voltage)) / length(voltage) * 1000 / dt
get_spike_rate(voltage::Matrix{T}) where {T<:Real} = get_spike_rate(voltage[1, :])
get_isi(spiketimes = Vector{Float32}) = diff(spiketimes)
cap_voltage(voltage::Vector{T}) where {T<:Real} =
    [v > postspike.AP_membrane ? -50 : v for v in voltage]

function get_cv(voltage)
    spiketimes = get_spike_times(voltage)
    intervals = get_isi(spiketimes)
    _cv = sqrt(var(intervals) / mean(intervals)^2)
    return isnan(_cv) ? 0.0 : _cv
end

function CV_isi(intervals::Vector{Float32})
    # input::Union{Matrix{Float32}, Vector{Float32}, BitArray{1}})
    _cv = sqrt(var(intervals) / mean(intervals)^2)
    return isnan(_cv) ? 0.0 : _cv
end


function get_cv2(voltage)
    spiketimes = get_spike_times(voltage)
    intervals = get_isi(spiketimes)

    _cv = sqrt(var(intervals) / mean(intervals)^2)
    return isnan(_cv) ? 0.0 : _cv
end
