using Main
module TripodNeuron
using Logging, Printf, Distributions
# include("TripodNeuron/base/structs/dendrites.jl")
include("../TripodNeuron/base.jl")
postspike = PostSpike(10, round(Int, 1 / dt), 1 / 50)
AdEx = get_bonoclopath_params()
end

TN = TripodNeuron
using Plots
using LinearAlgebra
using StatsBase


include(joinpath(@__DIR__, "STDP_rules.jl"))

stdp = bono_vstdp()
##


function protocol(;
    tripod = nothing,
    step::Int64,
    dt::Float32,
    AMPAsynapses::Matrix{Float32},
    NMDAsynapses::Matrix{Float32},
    spikes::BitArray{3},
    inh_rate::Float32,
)

    for x in findall(spikes[:, step, 1])
        TN.exc_spike_plastic!(
            tripod.d[1],
            AMPA_eff = AMPAsynapses[x, 1],
            NMDA_eff = NMDAsynapses[x, 1],
        )
    end
    for x in findall(spikes[:, step, 2])
        TN.exc_spike_plastic!(
            tripod.d[2],
            AMPA_eff = AMPAsynapses[x, 2],
            NMDA_eff = NMDAsynapses[x, 2],
        )
    end
    TN.set_rate!(tripod.d[1], -inh_rate)
    TN.set_rate!(tripod.d[2], -inh_rate)
    TN.set_rate!(tripod.s, 1.0, eff = 1.5)
    TN.set_rate!(tripod.s, -0.5, eff = 1.0)
end


function simulate_learning(
    t::TN.Tripod,
    simTime::Int64,
    stimulation_protocol::Union{Function,Array{Int64}};
    NMDAsynapses::Matrix{Float32},
    AMPAsynapses::Matrix{Float32},
    spikes,
    kwargs...,
)
    total_steps = round(Int, simTime / TN.dt)
    n_dend = length(t.d)
    ## Recordings
    currents = Array{Float32,1}(undef, n_dend + 1)
    voltage = Array{Float32,2}(undef, n_dend + 1, total_steps)
    spiked = false
    last_spike = -100.0
    spiked = (t.s.v > TN.AdEx.θ)

    θ = Array{Float32,1}(undef, total_steps)
    u⁻ = zeros(2)
    u⁺ = zeros(2)

    xstdp = zeros(size(AMPAsynapses))
    delay_v = zeros(2, round(Int, 1 / TN.dt))
    trace_xstdp = zeros(size(AMPAsynapses)..., total_steps)
    rec_synapses = zeros(size(AMPAsynapses)..., total_steps)
    initial_synapses = sum(AMPAsynapses, dims = 1)[1, :]
    norm = zeros(2)

    @fastmath @inbounds for tt = 1:total_steps


        if TN.update_tripod!(t, currents, spiked)
            spiked = true
            last_spike = tt
        else
            spiked = false
        end
        xstdp .+= spikes[:, tt, :]
        for d in [1, 2]
            delay_v[d, 1:end-1] = delay_v[d, 2:end]
            delay_v[d, end] = t.d[d].v
        end


        #############################
        ## Store results
        voltage[1, tt] = t.s.v
        for n = 1:n_dend
            voltage[1+n, tt] = t.d[n].v
        end
        θ[tt] = t.s.θ
        #############################

        stimulation_protocol(;
            tripod = t,
            step = tt,
            dt = TN.dt,
            AMPAsynapses = AMPAsynapses,
            NMDAsynapses = NMDAsynapses,
            spikes = prespikes,
            kwargs...,
        )

        if (tt % 200 == 0)
            norm[:] .= initial_synapses ./ (sum(AMPAsynapses, dims = 1)[1, :])
            rmul!((@view AMPAsynapses[:, 1]), norm[1])
            rmul!((@view AMPAsynapses[:, 2]), norm[2])
        end

        # @TODO in the BonoClopath there is an additional time-scale τ₁.

        ## STDP Plasticity Rules
        xstdp .+= -TN.dt .* (xstdp) .* stdp.τx⁻    #x

        # synaptic plasticity depends on 'pre' synaptic neurons
        #LTD
        for d in [1, 2]
            u⁻[d] += TN.dt * (delay_v[d, 1] - u⁻[d]) * stdp.τu⁻ #u
            u⁺[d] += TN.dt * (delay_v[d, 1] - u⁺[d]) * stdp.τv⁻ #v
            #LTD
            if u⁻[d] > stdp.θ⁻
                AMPAsynapses[findall(prespikes[:, tt, d]), d] .-=
                    TN.dt * stdp.a⁻ * (u⁻[d] - stdp.θ⁻)
            end
            #LTP
            if (u⁺[d] > stdp.θ⁻) && (t.d[d].v > stdp.θ⁺)
                AMPAsynapses[:, d] .+=
                    TN.dt .* xstdp[:, d] .* stdp.a⁻ *
                    (t.d[d].v - stdp.θ⁺) *
                    (u⁺[d] - stdp.θ⁻)
            end
        end
        AMPAsynapses[findall(AMPAsynapses .> stdp.j⁺)] .= stdp.j⁺
        AMPAsynapses[findall(AMPAsynapses .< stdp.j⁻)] .= stdp.j⁻
        rec_synapses[:, :, tt] .= AMPAsynapses
        trace_xstdp[:, :, tt] .= xstdp
    end
    ## Current and Voltage plot
    return voltage, rec_synapses, trace_xstdp
end



function make_spikes(
    simtime::Int64,
    dt::Float32;
    n_synapses::Int64,
    distal_rate::Float32,
    proximal_rate::Float32,
    interval_duration::Int64,
    swap_order::Bool,
)

    steps = round(Int, simtime / dt)
    interval_steps = round(Int, interval_duration / dt)
    intervals = floor(Int, simtime / interval_duration)

    pre_neurons = collect(1:n_synapses)

    distal = falses(n_synapses, steps)
    proximal = falses(n_synapses, steps)
    states = zeros(steps)
    for t = 1:intervals
        interval = 1+interval_steps*(t-1):interval_steps*t
        states[interval] .= 0
        # for s in 1:n_synapses-1
        # 	distal[s, interval] = TN.PoissonInput(proximal_rate/n_synapses,interval_duration,TN.dt)
        # end
        if t % 6 == 1
            states[interval] .= 1
            proximal[1:3, interval] =
                TN.PoissonInput(proximal_rate, interval_duration, TN.dt, neurons = 3)
            if !swap_order
                distal[:, interval] = TN.PoissonInput(
                    distal_rate,
                    interval_duration,
                    TN.dt,
                    neurons = n_synapses,
                )
            end
        elseif t % 6 == 3
            states[interval] .= 3
            proximal[4:6, interval] =
                TN.PoissonInput(proximal_rate, interval_duration, TN.dt, neurons = 3)
            if swap_order
                distal[:, interval] = TN.PoissonInput(
                    distal_rate,
                    interval_duration,
                    TN.dt,
                    neurons = n_synapses,
                )
            end
        elseif t % 6 == 5
            states[interval] .= 5
            proximal[7:9, interval] =
                TN.PoissonInput(proximal_rate, interval_duration, TN.dt, neurons = 3)
        end
    end
    return distal, proximal, states
end
##
using EllipsisNotation
using RollingFunctions

SIMTIME = 5000
IN_NEURONS = 10
PROX_RATE = 800.0
DIST_RATE = 1500.0
INTERVAL = 80
REPETITION = 2
EFF = stdp.j⁻ * 1.5
INH_RATE = 2.0

t = TN.Tripod(TN.H_distal_proximal)
_ = TN.simulate_fast(t, TN.null_input, 1000)
distal, proximal, states = make_spikes(
    SIMTIME,
    TN.dt,
    n_synapses = IN_NEURONS,
    distal_rate = DIST_RATE,
    proximal_rate = PROX_RATE,
    interval_duration = INTERVAL,
    swap_order = false,
)
prespikes = cat(distal, proximal, dims = 3)
AMPAsynapses_start = EFF .* ones(IN_NEURONS, 2)
NMDAsynapses = EFF .* ones(IN_NEURONS, 2)
voltage, AMPAsynapses, trace = simulate_learning(
    t,
    SIMTIME,
    protocol,
    AMPAsynapses = AMPAsynapses_start,
    NMDAsynapses = NMDAsynapses,
    spikes = prespikes,
    inh_rate = INH_RATE,
)

distal, proximal, states_2 = make_spikes(
    SIMTIME,
    TN.dt,
    n_synapses = IN_NEURONS,
    distal_rate = DIST_RATE,
    proximal_rate = PROX_RATE,
    interval_duration = INTERVAL,
    swap_order = true,
)
prespikes = cat(distal, proximal, dims = 3)
presynapses = AMPAsynapses[:, :, end]

voltage_2, AMPAsynapses_2, trace_2 = simulate_learning(
    t,
    SIMTIME,
    protocol,
    AMPAsynapses = AMPAsynapses_start,
    NMDAsynapses = NMDAsynapses,
    spikes = prespikes,
    inh_rate = INH_RATE,
)
voltage = cat(voltage, voltage_2, dims = 2)
AMPAsynapses = cat(AMPAsynapses, AMPAsynapses_2, dims = 3)
states = cat(states, states_2, dims = 3)

##
p = plot(runmean(mean(AMPAsynapses[1:3, 2, :], dims = 1)[1, :], 700), label = "Syn 1")
p = plot!(runmean(mean(AMPAsynapses[4:6, 2, :], dims = 1)[1, :], 700), label = "Syn 2")
p = plot!(
    runmean(mean(AMPAsynapses[7:9, 2, :], dims = 1)[1, :], 700),
    label = "Syn 3",
    ylabel = "Synaptic \nEfficacy",
)
p = plot!(xticks = (0:25000:100000, 0:2.5:10))
# p = plot!(mean(AMPAsynapses[1:3,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[4:6,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[7:9,1,:], dims=1)[1,:], ls=:dash)
plot(4 * states[1:10000], c = :black, lw = 3, label = "")
plot!(zeros(10000), c = :white, lw = 4, label = "")
q = plot!(
    ylims = (-80, 40),
    voltage[:, 1:10000]',
    c = [:black :red :blue],
    label = false,
    ylabel = "Membrane \nPotential(mV)",
    xticks = (0:5000:10000, [0, 0.5, 1.0]),
    yticks = (-80:20:0, -80:20:0),
)
plot(4 * states[end-10000:end], c = :black, lw = 3, label = "")
plot!(zeros(10000), c = :white, lw = 4, label = "")
s = plot!(
    ylims = (-80, 40),
    voltage[:, end-10000:end]',
    c = [:black :red :blue],
    label = false,
    xticks = (0:5000:10000, [9.0, 9.5, 10]),
    yaxis = false,
)
layout = @layout [
    a{0.5h}
    a{0.45w} c{0.1w} a{0.45w}
]
plot(p, q, plot(frame = :none, xlabel = "\n Time (s)"), s, layout = layout)
