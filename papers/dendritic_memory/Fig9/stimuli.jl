
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using HDF5, Random, RollingFunctions, Distributions, StatsBase


"""
Oscillatory signal with ´interval´ (ms) duration (half wave-length). 
It starts in state ´0´

# Arguments
- `simtime` : the duration of the simulation
- `interval` : the duration of the oscillation
"""
function oscillating_signal(simtime, interval, dt=TN.dt)
    sig = zeros(Int,round(Int,simtime/dt))
    for step in 1:length(sig)
        fase =(step ÷ (interval/dt))
        if fase%2 == 0
                sig[step] = 0
        elseif fase%2 ==1
                sig[step] = 1
        end
    end
    return sig
end


"""
Produce signal with irregular interval duration. 
Interval length is drawn from an exponential distribution with variable interval length. 
It starts in state ´0´

# Arguments
- `simtime` : the duration of the simulation
- `interval` : the average duration of the interval
"""
function exp_decay_signal(simtime, interval, dt=TN.dt)
    signal = zeros(round(Int,simtime/dt))
    mem    = 0
    next   = 1
    for tt in 1:simtime*10
        if tt == next
            next = round(Int,rand(Exponential(interval/dt))+ tt+1)
            mem = (mem+1)%2
        end
        signal[tt] = mem
    end
    signal = sign.(signal)
    return signal
end

"""
Three-phasic signal function with ´interval´ (ms) duration (half wave-length). 
It starts in state ´0´
# Arguments
- `simtime` : the duration of the simulation
- `interval` : the duration of the oscillation
"""
function three_phase_signal(simtime, interval, dt=TN.dt)
    sig = zeros(Int,round(Int,simtime/dt))
    for step in 1:length(sig)
        fase =(step ÷ (interval/dt))
        if fase%3 == 0
                sig[step] = 0
        elseif fase%3 ==1
                sig[step] = 1
        elseif fase%3 ==2
                sig[step] = 2
        end
    end
    return sig
end

"""
Extract the time intervals where the stimulus is in 
"""
function get_stimuli(signal; states=2, dt=0.1)
    traces = []
    for n in 1:states
        push!(traces, findall(signal .==n-1))
    end
    return traces
end

function get_stimuli(signal; states=2, dt=0.1)
    traces = []
    for n in 1:states
        push!(traces, findall(signal .==n-1))
    end
    return traces
end
"""
Apply stimuli to the input spike train
"""
function apply_stimuli!(stimuli, inputs; soma=true)
    inputs[1,:] = TN.PoissonInput(500,length(inputs[1,:]), 1.f0)
    ## for stim in 1: increase excitation on dendrite A
    inputs[[2],stimuli[1]] .*=2.
    ## for stim in 2: increase inhibition on dendrite B
    inputs[[3],stimuli[2]] .*=2.
end

"""
Apply stimuli to the input spike train with non balanced condition
"""
function apply_three_phase_stimuli!(stimuli, inputs, order)
    ## Very high classification score with this input configuration
    # inputs[:] .*= 2
    # inputs[:] .= 0
    # if order == 1
    #     inputs[2,stimuli[1]] = TN.PoissonInput(500,length(stimuli[1]), 1.f0)
    #     inputs[3,stimuli[2]] = TN.PoissonInput(100,length(stimuli[2]), 1.f0)
    # elseif order == 2
    #     inputs[3,stimuli[1]] = TN.PoissonInput(500,length(stimuli[1]), 1.f0)
    #     inputs[2,stimuli[2]] = TN.PoissonInput(100,length(stimuli[2]), 1.f0)
    # end
    ## Lower classification score, but consistent with balanced inputs
    inputs[1,:] = TN.PoissonInput(100,length(inputs[1,:]), 1.f0)
    if order == 1
        inputs[2,stimuli[1]] .*= 2.
        inputs[3,stimuli[2]] .*= 2.
    elseif order == 2
        inputs[2,stimuli[2]] .*= 2.
        inputs[3,stimuli[1]] .*= 2.
    end
    inputs[2,stimuli[3]] .*= 0.3
    inputs[3,stimuli[3]] .*= 0.3
end

function get_label(simtime, stimuli; order)
    labels = zeros(Int, round(Int,simtime/TN.dt))
    if order == 1
        labels[stimuli[1]] .= 1 
        labels[stimuli[2]] .= 2
        labels[stimuli[3]] .= 3
    else
        labels[stimuli[1]] .= 2 
        labels[stimuli[2]] .= 1
        labels[stimuli[3]] .= 3
    end
    return labels
end


#==============================================================================
            Measure signal/spike rate correlations
==============================================================================#

@inline function exp32(x::Real)
    x = ifelse(x < -10f0, -32f0, x)
    x = 1f0 + x / 32f0
    x *= x; x *= x; x *= x; x *= x; x *= x
    return x
end

Θ(x::Float64) = x > 0.0 ? x : 0.0

function alpha_function(t; t0, τ)
	if abs(t-t0)/τ > 5
		return 0.f0
	else
    return (t - t0) / τ * exp32(1 - (t - t0) / τ) * Θ(1.0 * (t - t0))
	end
end

"""
Convolve with alpha function of timescale τ
"""
function convolve(spiketime::Vector; interval::AbstractRange, τ = 100)
    rate = zeros(Float32,length(interval))
    for i in eachindex(interval)
        v = 0
         @simd for t0 in spiketime
         @fastmath v += alpha_function(interval[i], t0 = t0, τ = τ)
        end
        rate[i] = v/τ
    end
    return rate
end


""""
Integrate the spike train into a rate by convolving each spike with an alpha function.
"""
function spikes_integrator(model, signal, simtime; samples=200)
    spikes = Vector{}()
    stimuli = get_stimuli(signal)
    spiketimes = Vector{Vector{Float32}}()
    for _ in 1:samples
        inputs = TN.make_balance(simtime, model=model)
        apply_stimuli!(stimuli, inputs)
        v = TN.run_tripod(inputs, simtime; model...)
        append!(spikes,TN.get_spike_times(v[1,:]) *TN.dt)
        push!(spiketimes, TN.get_spike_times(v[1,:])* TN.dt)
    end
    return convolve(Float32.(spikes), interval=1:1.: simtime, τ=10), spiketimes
end


function signal_integrator(signal; interval, τ=10)
    events = findall(abs.(diff(signal[1:end])) .>0) 
    return convolve(events * 0.1, interval=interval, τ=τ)
end


function normalize(s)
    train_std = StatsBase.fit(ZScoreTransform, s[1000:end], dims=1)
    StatsBase.transform!(train_std,s[1000:end])
end

function measure_correlation(signal, output; max_shift=100)
    signal = normalize(signal)
    output = normalize(output)
    shifts = -max_shift+1:5:max_shift
    zz = zeros(length(shifts))
    for (n,s) in enumerate(shifts)
        xx0 = (max_shift+1:length(signal) - max_shift) .+s
        xx1 = max_shift+1:length(signal) - max_shift
        zz[n] = cor(signal[xx0], output[xx1])
    end
    return zz
end


#==============================================================================
            Run simulation with sequence structured input
==============================================================================#

function simulate_sequence(;model, order, β::Real=1, interval::Real=300, τinh::Real=50, repeat=5, A=2)
    simtime = minimum([interval*repeat*3, 10_000])
    spikes =  TN.make_balance(simtime, model= model)
    signal = three_phase_signal(simtime, interval)
    stimuli = get_stimuli(signal,states=3)
    apply_three_phase_stimuli!(stimuli, spikes, order)
    if sum(model.ds) == 0
        (spikes.*=3)
    end
    voltage, _ = TN.run_tripod_sequence(spikes, simtime, β=β, τinh=τinh, A=A; model...)
    labels = get_label(simtime, stimuli, order=order)
    return voltage[1,:], labels
end


function run_batch(; model, β, interval, τinh, batch_size=200, repeat =10, A=5 )
    rates = zeros(batch_size)
    cvs =    zeros(batch_size)
    orders = zeros(batch_size)
    Threads.@threads   for x in 1:batch_size
        order = rand([1,2])
        v, _ = simulate_sequence(;model=model, order=order, β=β, interval=interval, τinh=τinh, repeat=repeat, A=A)
        orders[x] = order
        cvs[x] = TN.get_cv(v)
        rates[x] = TN.get_spike_rate(v)
    end
    X = hcat(cvs, rates)
    Y = orders 
    return X, Y
end
