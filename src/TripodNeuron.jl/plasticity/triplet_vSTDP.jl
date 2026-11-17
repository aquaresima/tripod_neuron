using Main
using Plots
using EllipsisNotation

path_spiketimit = "/home/cocconat/Documents/Research/phd_project/speech/Spike TIMIT/src/SpikeTimit.jl"
include(path_spiketimit)
include("learning_functions.jl")
triplet = pfister_visualcortex()
vstdp = clopath_vstdp_visualcortex()

##
SIMTIME = 1000
IN_NEURONS = 1
PROX_RATE = 2.0


ρ = 36.5 # input rate
Δt1 = 10
Δt2 = 15
REPETITIONS = 5
EFF = 1.0
pre, post = make_triplet(SIMTIME, TN.dt, Δt1, Δt2, ρ)
SpikeTimit.matrix2sequence(pre, 0.1)

##

EFF = 1.0
neurons = 1
AMPAsynapses = EFF .* ones(neurons)
NMDAsynapses = EFF .* ones(neurons)
GABAsynapses = 0.0 .* ones(neurons)
pre_spikes = TN.PoissonInput(100, SIMTIME, TN.dt, neurons = neurons)
post_spikes = TN.PoissonInput(50, SIMTIME, TN.dt, neurons = 1)
all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre_spikes, TN.dt, "ms")
post_firingtimes, _ = SpikeTimit.matrix2sequence(post_spikes, TN.dt, "ms")

s = TN.AdEx_Soma()
w, v, θ = simulate_clopath(
    s,
    SIMTIME,
    AMPAsynapses,
    NMDAsynapses,
    GABAsynapses,
    all_ft_dt,
    all_neurons,
    post_firingtimes = post_firingtimes,
)

q = plot(v, ylabel = "Mem \nvoltage (mV)")
p = plot(θ, ylabel = "Firing
  \nthreshold (mV)")
s = plot(w[:, :]', ylabel = "Weight (pA)\n")
# scatter!(TN.get_spike_times(pre_spikes)[1],ones(51)*-40, ms=10)
# scatter!(TN.get_spike_times(post_spikes)[1],ones(51)*-40, ms=10)
# w
plot(p, s, q, layout = (3, 1), legend = false)
##
SIMITIME = 10000
IN_NEURONS = 1
plots = []
include("learning_functions.jl")
for EFF in [1.0, 20.0, 50.0, 100.0]
    GABAsynapses = zeros(0)
    NMDAsynapses = EFF .* ones(IN_NEURONS)
    pre = TN.PoissonInput(10, SIMTIME, TN.dt)
    post = TN.PoissonInput(1, SIMTIME, TN.dt)
    current = 1600.0

    AMPAsynapses = EFF .* ones(IN_NEURONS)
    all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
    post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
    s = TN.AdEx_Soma()

    w, v, θ = simulate_pfister(
        s,
        SIMTIME,
        AMPAsynapses,
        NMDAsynapses,
        GABAsynapses,
        all_ft_dt,
        all_neurons,
        post_firingtimes = post_firingtimes,
        max_efficacy = 1000.0,
        current = current,
    )

    plot(w[1, 1:50:end] ./ w[1, 1], label = "Triplet")

    AMPAsynapses = EFF .* ones(IN_NEURONS)
    all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
    post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
    s = TN.AdEx_Soma()

    w, v, θ = simulate_clopath(
        s,
        SIMTIME,
        AMPAsynapses,
        NMDAsynapses,
        GABAsynapses,
        all_ft_dt,
        all_neurons,
        post_firingtimes = post_firingtimes,
        max_efficacy = 1000.0,
        current = current,
    )

    p = plot!(
        w[1, 1:50:end] ./ w[1, 1],
        title = "Syn eff: " * string(EFF),
        labels = "vSTDP",
    )
    # plot!( xticks=(0:500:2000, 0:2.5:10))
    if EFF == 1.0
        plot!(legend = :topleft)
    else
        plot!(legend = false)
    end
    push!(plots, p)
end
plot(plots..., ylabel = "relative ΔW", xlabel = "Time(s)")

##
pre = TN.PoissonInput(10, SIMTIME, TN.dt)
post = TN.PoissonInput(1, SIMTIME, TN.dt)
for current in [0.0, 100.0, 500.0, 1000.0]
    EFF = 1
    GABAsynapses = zeros(0)
    NMDAsynapses = EFF .* ones(IN_NEURONS)

    AMPAsynapses = EFF .* ones(IN_NEURONS)
    all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
    post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
    s = TN.AdEx_Soma()

    w, v, θ = simulate_pfister(
        s,
        SIMTIME,
        AMPAsynapses,
        NMDAsynapses,
        GABAsynapses,
        all_ft_dt,
        all_neurons,
        post_firingtimes = post_firingtimes,
        max_efficacy = 1000.0,
        current = current,
    )

    plot(w[1, 1:50:end] ./ w[1, 1], label = "Triplet")

    AMPAsynapses = EFF .* ones(IN_NEURONS)
    all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
    post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
    s = TN.AdEx_Soma()

    w, v, θ = simulate_clopath(
        s,
        SIMTIME,
        AMPAsynapses,
        NMDAsynapses,
        GABAsynapses,
        all_ft_dt,
        all_neurons,
        post_firingtimes = post_firingtimes,
        max_efficacy = 1000.0,
        current = current,
    )

    p = plot!(
        w[1, 1:50:end] ./ w[1, 1],
        title = "current: " * string(current),
        labels = "vSTDP",
    )
    # plot!( xticks=(0:500:2000, 0:2.5:10))
    if EFF == 1.0
        plot!(legend = :topleft)
    else
        plot!(legend = false)
    end
    push!(plots, p)
end
plot(plots..., ylabel = "relative ΔW", xlabel = "Time(s)")
# plots[1]
##
pre_rates = 1:5:51
post_rates = 1:5:51
SIMTIME = 10000

EFF = 10.0
# function estimate_learning(pre_rates, post_rates)
Δw = zeros(SIMTIME, size(pre_rates)..., size(post_rates)..., 2)
Threads.@threads for n = 1:length(pre_rates)
    for m = 1:length(post_rates)
        NMDAsynapses = ones(IN_NEURONS)
        pr = pre_rates[n]
        po = post_rates[m]
        @show po, n
        @show pr, m
        pre = TN.PoissonInput(pr, SIMTIME, TN.dt)
        post = TN.PoissonInput(po, SIMTIME, TN.dt)

        AMPAsynapses = EFF .* ones(IN_NEURONS)
        all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
        post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
        s = TN.AdEx_Soma()
        w, v, θ = simulate_clopath(
            s,
            SIMTIME,
            AMPAsynapses,
            NMDAsynapses,
            GABAsynapses,
            all_ft_dt,
            all_neurons,
            post_firingtimes = post_firingtimes,
            max_efficacy = 10.0,
        )
        Δw[:, n, m, 1] .= w[1:10:end]

        AMPAsynapses = EFF .* ones(IN_NEURONS)
        all_ft_dt, all_neurons = SpikeTimit.matrix2sequence(pre, TN.dt, "ms")
        post_firingtimes, _ = SpikeTimit.matrix2sequence(post, TN.dt, "ms")
        s = TN.AdEx_Soma()
        w, v, θ = simulate_clopath(
            s,
            SIMTIME,
            AMPAsynapses,
            NMDAsynapses,
            GABAsynapses,
            all_ft_dt,
            all_neurons,
            post_firingtimes = post_firingtimes,
            max_efficacy = 10.0,
        )
        Δw[:, n, m, 2] .= w[1:10:end]
    end
end

Δw[end, ..]
##
# 	return weights
# end
# Δw = estimate_learning(pre_rates, post_rates)
values = mean(Δw, dims = 1)[1, :, :, :]
values = (Δw[end, :, :, :] .- 1)
plot(
    heatmap(
        pre_rates,
        post_rates,
        values[.., 1]',
        xlabel = "Pre rate (Hz)",
        ylabel = "Post rate (Hz)",
        title = "Triplet",
        c = :redsblues,
    ),
    heatmap(
        pre_rates,
        post_rates,
        values[.., 2]',
        clims = (-1, 1),
        xlabel = "Pre rate (Hz)",
        ylabel = "Post rate (Hz)",
        c = :redsblues,
    ),
    title = "vSTDP",
    clims = (-1, 1),
    layout = (2, 1),
)



##


function make_triplet(simtime, dt, Δt1, Δt2, ρ)
    prespikes = falses(1, round(Int, simtime / dt))
    postspikes = falses(round(Int, simtime / dt))
    interval = round(Int, 1 / ρ / dt * 1000)
    repetitions = round(Int, simtime / dt / interval)
    @show interval, repetitions
    for x = 1:repetitions-2
        prespikes[x*interval] = true
        prespikes[x*interval+round(Int, Δt1 / dt)] = true
        postspikes[x*interval+round(Int, (Δt1 + Δt2) / dt)] = true
    end
    return prespikes, postspikes
end
using EllipsisNotation
SIMTIME = 10000
IN_NEURONS = 1
PROX_RATE = 2.0
ρ = 36.5 # input rate
Δt1 = 10
Δt2 = 15
REPETITIONS = 5
EFF = 1.0
#
# EFF = 1.5
plots = []
for EFF in [1.0, 20.0, 50.0, 100.0]
    AMPAsynapses_start = EFF .* ones(IN_NEURONS)
    NMDAsynapses = EFF .* ones(IN_NEURONS)

    pre = TN.PoissonInput(20, SIMTIME, TN.dt)
    post = TN.PoissonInput(50, SIMTIME, TN.dt)[1, :]
    s = TN.AdEx_Soma()
    w = simulate_learning_pfister(
        s,
        SIMTIME,
        protocol_triplet,
        AMPAsynapses = AMPAsynapses_start,
        NMDAsynapses = NMDAsynapses,
        prespikes = pre,
        postspikes = post,
    )
    plot(w[1:50:end] ./ w[1], label = "Triplet")
    s = TN.AdEx_Soma()
    w = simulate_learning_clopath(
        s,
        SIMTIME,
        protocol_triplet,
        AMPAsynapses = AMPAsynapses_start,
        NMDAsynapses = NMDAsynapses,
        prespikes = pre,
        postspikes = post,
        max_efficacy = 1000.0,
    )
    # mean(w)
    p = plot!(
        w[1:50:end] ./ w[1],
        title = "Syn eff: " * string(EFF),
        labels = "vSTDP",
        xticks = (0:500:2000, 0:2.5:10),
    )
    push!(plots, p)
end

using LaTeXStrings
plot(plots..., ylabel = "Relative " * L"\Delta" * "w", xlabel = "Time (s)")


##
pre_rates = 1:2:41
post_rates = 1:2:41
##
# function estimate_learning(pre_rates, post_rates)
Δw = zeros(1000, size(pre_rates)..., size(post_rates)...)
for n = 1:length(pre_rates)
    for m = 1:length(post_rates)
        s = TN.AdEx_Soma()
        AMPAsynapses_start = ones(IN_NEURONS)
        NMDAsynapses = ones(IN_NEURONS)
        pr = pre_rates[n]
        po = post_rates[m]
        @show po, n
        @show pr, m
        pre = TN.PoissonInput(pr, SIMTIME, TN.dt)
        post = TN.PoissonInput(po, SIMTIME, TN.dt)[1, :]
        w = simulate_learning_pfister(
            s,
            SIMTIME,
            protocol_triplet,
            AMPAsynapses = AMPAsynapses_start,
            NMDAsynapses = NMDAsynapses,
            prespikes = pre,
            postspikes = post,
        )
        Δw[:, n, m] .= w[1:100:end]
    end
end
# 	return weights
# end
# Δw = estimate_learning(pre_rates, post_rates)
heatmap(
    pre_rates,
    post_rates,
    mean(Δw, dims = 1)[1, :, :]' .- 1,
    clims = (-1, 1),
    xlabel = "Pre rate (Hz)",
    ylabel = "Post rate (Hz)",
    c = :redsblues,
)
##
n, m = 1, 3
s = TN.AdEx_Soma()
AMPAsynapses_start = ones(IN_NEURONS)
NMDAsynapses = ones(IN_NEURONS)
pr = pre_rates[n]
po = post_rates[m]
@show po, pr
pre = TN.PoissonInput(pr, SIMTIME, TN.dt)
post = TN.PoissonInput(po, SIMTIME, TN.dt)[1, :]
w = simulate_learning_soma(
    s,
    SIMTIME,
    protocol_triplet,
    AMPAsynapses = AMPAsynapses_start,
    NMDAsynapses = NMDAsynapses,
    prespikes = pre,
    postspikes = post,
)
plot(Δw[:, n, m])
plot!(w[1:50:end])
##

plot(voltage)
plot(w)

# plot(rec_triplet', label=["r1" "r2" "o1" "o2"])
##
p = plot(mean(AMPAsynapses[1:3, 2, :], dims = 1)[1, :])
p = plot!(mean(AMPAsynapses[4:6, 2, :], dims = 1)[1, :])
p = plot!(mean(AMPAsynapses[7:9, 2, :], dims = 1)[1, :])
# p = plot!(mean(AMPAsynapses[1:3,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[4:6,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[7:9,1,:], dims=1)[1,:], ls=:dash)
q = plot(voltage', c = [:black :red :blue])
plot(p, q, layout = (2, 1))
# plot(states)
##

distal, proximal, states = make_spikes(
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
    inh_rate = 1.0,
)
voltage = cat(voltage, voltage_2, dims = 2)
AMPAsynapses = cat(AMPAsynapses, AMPAsynapses_2, dims = 3)
##
using RollingFunctions

p = plot(runmean(mean(AMPAsynapses[1:3, 2, :], dims = 1)[1, :], 700))
p = plot!(runmean(mean(AMPAsynapses[4:6, 2, :], dims = 1)[1, :], 700))
p = plot!(runmean(mean(AMPAsynapses[7:9, 2, :], dims = 1)[1, :], 700))
# p = plot!(mean(AMPAsynapses[1:3,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[4:6,1,:], dims=1)[1,:], ls=:dash)
# p = plot!(mean(AMPAsynapses[7:9,1,:], dims=1)[1,:], ls=:dash)
q = plot(voltage[:, 1:15000]', c = [:black :red :blue])
s = plot(voltage[:, end-15000:end]', c = [:black :red :blue])
layout = @layout [
    a{0.3h}
    grid(1, 2)
]
plot(p, q, s, layout = layout)
##
# plot(states)
#
#  synapses
# full_volt = zeros(size(voltage)..., REPETITION*SIMTIME)
# full_syn = zeros(3, REPETITION*SIMTIME)
# # full_syn = zeros(size()..., 5)
# for i in 1:REPETITION
# 	distal, proximal, states2 = make_spikes(SIMTIME, TN.dt, n_synapses=IN_NEURONS,
#  	distal_rate=DIST_RATE, proximal_rate=PROX_RATE, interval_duration=INTERVAL, swap_order=false)
# 	prespikes = cat(distal, proximal, dims=3)
#
# 	voltage, synapses, trace= simulate_learning(t, SIMTIME, protocol, synapses =synapses[:,:,end,1], spikes=prespikes, inh_rate=2.)
# 	for i in 1:3
# 	full_syn[..,SIMTIME*(i-1)+1:SIMTIME*i] .= mean(full_syn[i-1*3+1:3,2,:,1], dims=1)[1,:])
# p = plot!(mean(full_syn[4:6,2,:,1], dims=1)[1,:])
# p = plot!(mean(full_syn[7:9,2,:,1], dims=1)[1,:])
# end
# #
# # full_syn = cat(synapses, synapses2, dims=3)
# # full_volt = cat(voltage, voltage2, dims=2)
# # full_states = cat(states[:], states2[:], dims=1)
# p  = plot(mean(full_syn[1:3,2,:,1], dims=1)[1,:])
# p = plot!(mean(full_syn[4:6,2,:,1], dims=1)[1,:])
# p = plot!(mean(full_syn[7:9,2,:,1], dims=1)[1,:])

##
plotlyjs()
p = plot(synapses[:, 2, :, 1]', labels = "")
# p = plot!(synapses[:,1,:,2]', labels="")
plot!(trace[:, 1, :]' ./ maximum(trace), ls = :dot, labels = "")
# scatter(distal[:,:]')
q = plot(voltage', c = [:black])
plot(p, q, layout = (2, 1))
##
## plot(full_volt ./80 .+ 1.5, c=:black)
# plotlyjs()
##
# q = plot(full_syn[1,2,:,1])
# plot!(full_syn[2,2,:,1])
# plot!(full_syn[3,2,:,1])

p = plot(
    full_volt[1, :],
    xlims = (1, 8000),
    c = :black,
    xticks = (0:2000:8000, string.(0:200:800)),
    labels = "",
)
p = plot(
    full_volt[1, :],
    xlims = (96000 - 7999, 96000),
    c = :black,
    xticks = (96000-8000:2000:96000, string.(0:200:800)),
    labels = "",
)
p = plot(full_volt[1, :], xlims = (1, 8000), c = :black)
plot(p, q, p, labels = "", layout = [a{0.3w}])
##
# pyplot()
plot(full_volt[1, :], c = :black)
plot!(full_volt[2, :])
plot!(full_volt[3, :])
plot!(full_states)

# plot(sum(full_syn[:,:,:,1], dims=[1,2])[1,1,:])
##

# plot!(mean(full_syn[:,1,:], dims=1)[1,:])

##
##
plotlyjs()
t = TN.Tripod(TN.H_distal_proximal)
IN_NEURONS = presynapses = 1.5 .* ones(IN_NEURONS, 2)
prespikes = cat(
    TN.PoissonInput(100, 1000, TN.dt, neurons = IN_NEURONS),
    TN.PoissonInput(100, 1000, TN.dt, neurons = IN_NEURONS),
    dims = 3,
)
voltage, θ, synapses =
    simulate_learning(t, 1000, protocol, synapses = presynapses, spikes = prespikes)

plot(plot(synapses[:, 1, :]'), plot(synapses[:, 2, :]'), layout = (2, 1))

plot(voltage[1, :], c = :black)
plot!(voltage[2, :], c = :red)
plot!(voltage[3, :], c = :green)


# plot([runmean(sum(proximal, dims=1)[1,:],1000)])
# plot!([runmean(sum(distal, dims=1)[1,:],1000)])
#  # sum(distal, dims=1)'])
# TN.raster_plot(TN.get_spike_times(proximal))
# TN.raster_plot(TN.get_spike_times(distal))
#
#
# TN.PoissonInput(100.,50,TN.dt, neurons=10)
