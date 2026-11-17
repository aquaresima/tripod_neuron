# Load Tripod
cd(@__DIR__)
using Main
module TripodNeuron
using Logging, Printf, Distributions
# include("TripodNeuron/base/structs/dendrites.jl")
include("../../TripodNeuron/base.jl")
postspike = PostSpike(10, round(Int, 1 / dt), 1 / 50)
AdEx = get_bonoclopath_params()
get_synapses_exc(dt, comp) = eyal_exc_synapse(dt, comp)
Esyn_dend = get_synapses_exc(dt, "dend")
Esyn_soma = get_synapses_exc(dt, "soma")
Isyn_sst = get_synapses_sst(dt)
Isyn_pv = get_synapses_pv(dt)
end


include("../learning_functions.jl")
include("../get_spiketimit.jl")
##
TN = TripodNeuron
using Plots
using LinearAlgebra
using StatsBase
import StatsBase.var
using RollingFunctions
using Distributions

dt = TN.dt

##
β = 0.03
function make_rates(sim_time::Int, dt::Float32; baserate::Int = 5, peakrate::Int = 500)
    """
    Produce inpute rates, dimension is kHz⁻¹
    Resolution is ms.
    """
    remove_second = true
    backgroundNoise = baserate / 1000
    noise = 0.0
    inputs = zeros(round(Int, sim_time))
    for t = 1:round(Int, sim_time)
        re = rand() - 0.5
        noise = re - (re - noise) * exp(-dt / 50)
        inputs[t] = backgroundNoise + maximum([0, noise * peakrate / 1000]) / β
    end
    return inputs
end


plot(make_rates(1000, 0.1) .* 100, ylabel = "Excitatory inputs (Hz)", xlabel = "Time (ms)")

function make_inputs(rate::Vector{Float32}, neurons::Int, dt::Float32)
    """
    Get input rates in kHz⁻¹ and produce 'neurons' spike-trains.
    Spike trains have resolution in 0.1 ms
    """
    # spike_train = falses(neurons, length(rate))
    spiketimes = Array{Vector{Float32}}(undef, neurons)
    timesteps = round(Int, length(rate) / dt)
    refractory = 0.0
    for n = 1:neurons
        spiketimes[n] = Vector{Float32}()
        for tt = 1:timesteps
            refractory -= dt
            if refractory < 0
                if rand() < rate[ceil(Int, tt * dt)] * dt
                    push!(spiketimes[n], tt * dt)
                    # refractory = 5
                end
            end
        end
    end
    return spiketimes
end


# SpikeTimit.raster_plot(exc_spiketimes).plt


##

function prepare_inputs(
    groups::Int,
    exc_neurons::Int,
    in_neurons::Int,
    sim_time::Int;
    baserate = 5,
    peakrate = 500,
)

    exc_spikes = SpikeTimit.Spiketimes()
    inh_spikes = SpikeTimit.Spiketimes()

    for x = 1:groups
        rate = make_rates(sim_time, 0.1, baserate = baserate, peakrate = peakrate)

        spikes_ = make_inputs(rate, exc_neurons, dt)
        exc_spikes = vcat(exc_spikes, spikes_)
        spikes_ = make_inputs(rate, in_neurons, dt)
        inh_spikes = vcat(inh_spikes, spikes_)
    end
    return exc_spikes, inh_spikes
end

function prepare_synapses(p::Int, groups::Int, exc_neurons::Int, in_neurons::Int)
    AMPAsynapses = ones(exc_neurons * groups)
    GABAsynapses = ones(in_neurons * groups)
    for x = 1:8
        We = 0.3 + 1.1 / (1 + (x - p)^4) #+ (rand() - 0.5)
        Wi = 0.1
        AMPAsynapses[1+(x-1)*exc_neurons:exc_neurons*x] .= We
        GABAsynapses[1+(x-1)*in_neurons:in_neurons*x] .= Wi
    end
    return AMPAsynapses, GABAsynapses
end
##

EXC_NEURONS = 100
IN_NEURONS = 25
sim_time = 1_000 * 10
GROUPS = 8
P = 5

exc_spiketimes, inh_spiketimes =
    prepare_inputs(GROUPS, EXC_NEURONS, IN_NEURONS, sim_time; baserate = 5, peakrate = 500)
AMPAsynapses, GABAsynapses = prepare_synapses(P, GROUPS, EXC_NEURONS, IN_NEURONS)
NMDAsynapses = 1 * ones(length(exc_spiketimes))
# plot(AMPAsynapses)

# plot(GABAsynapses)

##


##
# plot(length.(exc_spiketimes)./sim_time)
# # exc_spiketimes
# p= SpikeTimit.raster_plot(exc_spiketimes).plt;
# plot!(p, xlims=(1,2))
# plot(AMPAsynapses)
##
include("../learning_functions.jl")
include("../get_spiketimit.jl")

all_ft, all_neurons = SpikeTimit.spiketimes2sequence(exc_spiketimes, inh_spiketimes)
all_ft_dt = SpikeTimit.all_ft_to_dt(all_ft, TN.dt, "ms")

s = TN.AdEx_Soma()
iSTDP = vogels_istdp()
weights, voltage, θ, firing_rate = simulate_vogels(
    s,
    sim_time,
    NMDAsynapses,
    AMPAsynapses,
    GABAsynapses,
    all_ft_dt,
    all_neurons,
)


# plot(firing_rate)
plot(voltage)
plot(runmean(firing_rate, 1000))
##
epoch = round(Int, sim_time / 50)
heatmap(weights[:, 1:epoch:end], ylabel = "Inhibitory neurons", xlabel = "Epochs (200 s)")

##
ampas = zeros(8)
gabas = zeros(8)
for x = 1:8
    ampas[x] = mean(AMPAsynapses[1+(x-1)*100:x*100])
    gabas[x] = mean(weights[1+(x-1)*25:x*25, end])
end

plot(ampas, xlabel = "Neuron group ", label = "AMPA")
plot!(gabas, label = "GABA")
