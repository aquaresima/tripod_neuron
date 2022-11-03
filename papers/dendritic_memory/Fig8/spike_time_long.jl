"""
Test the encoding-recall protocol in the distal-proximal compartment.
The protocol is valid against a broad set of parameters, such a balanced condition (see comment in protocol)

This file reproduces the figures in section Dendritic memory (Fig8D, Fig 8E)
"""



using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using StatsBase, Statistics, HDF5 
using MLDataUtils, MLJLinearModels, RollingFunctions
using ProgressBars, Logging
##


"""
Protocol for memory encoding and access.

# Arguments

- `distal_rate` : The spike rate of the encoding signal
- `inh_rate` : The spike rate of the inhibition signal
- `inh_comp` : Compartment targeted by inhibition
- `probe` : Interval between encoding and access
- `proximal_rate` : Access cue signal (1kHz)
- `encoding_duration` : Duration of encoding phase
- `base_soma_rate` : background noise for the soma

"""
function protocol_memory(simtime;
    distal_rate::Real, 
    inh_rate::Real, 
    inh_comp::String,
    probe::Real,
    proximal_rate::Real = 1k,
    encoding_duration=ENCODING_TIME::Real= 50, 
    base_soma_rate = BASE_SOMA_RATE,
    model )

    spikes = TN.make_spikes(simtime, rate=0)
    ## base input on the soma
    spikes[1,:] = TN.PoissonInput(base_soma_rate, simtime, TN.dt, neurons=1) 
    spikes[2,:] = TN.PoissonInput(distal_rate, simtime, TN.dt, neurons=1)   * 5
    if inh_comp == "soma"
        spikes[4,:] = TN.PoissonInput(inh_rate, simtime, TN.dt, neurons=1)  * 2
    elseif inh_comp == "dend"
        spikes[5,:] = TN.PoissonInput(inh_rate, simtime, TN.dt, neurons=1)  * 2
    end
    spikes[3,:] = TN.PoissonInput(proximal_rate, simtime, TN.dt, neurons=1) * 2
    encoding_dt = round(Int, encoding_duration/TN.dt) 
    probe_dt = round(Int, probe/TN.dt) 
    # # encoding phase
    spikes[3,1: encoding_dt+probe_dt] .=0
    # # recall phase
    spikes[2,encoding_dt+1:end] .=0
    spikes[4,encoding_dt+1:end] .=0
    spikes[5,encoding_dt+1:end] .=0
    return spikes
end

# """ 
# Protocol with excitation-inhibition balance.
# This protocol also can be used, the results are slightly different also show the presence of a slow decaying dendritic memory 
# """
# function protocol_memory_balanced(simtime;
#                 distal_rate::Real, 
#                 inh_rate::Real, 
#                 inh_comp::String,
#                 probe::Real,
#                 proximal_rate::Real = 1k,
#                 encoding_duration=ENCODING_TIME::Real= 50, 
#                 base_soma_rate = BASE_SOMA_RATE,
#                 model)

#     encoding_dt = round(Int, encoding_duration/TN.dt) 
#     probe_dt = round(Int, probe/TN.dt) 
#     spikes = TN.make_balance(simtime, model=model)
#     ## base input on the soma
#     spikes[1,:] = TN.PoissonInput(base_soma_rate, simtime, TN.dt, neurons=1) 
#     spikes[2,1:encoding_duration] = TN.PoissonInput(distal_rate,encoding_duration, 1.f0, neurons=1)   * 5
#     if inh_comp == "soma"
#         spikes[4,1:encoding_duration] = TN.PoissonInput(inh_rate,encoding_duration, 1.f0, neurons=1)  * 2
#     elseif inh_comp == "dend"
#         spikes[5,1:encoding_duration] = TN.PoissonInput(inh_rate, encoding_duration, 1.f0, neurons=1)  * 2
#     end
#     spikes[3,:] = TN.PoissonInput(proximal_rate, simtime, TN.dt, neurons=1) * 10
#     spikes[3,1: encoding_dt+probe_dt] .=0
#     return spikes
# end




"""
Measure first spike time (MSF) for different conditions:

# Arguments
- `simtime` : the duration of the simulation
- `distal_rate` : the spike rate of the encoding signal
- `inh_rate` : The spike rate of the inhibition signal
- `inh_comp` : Compartment targeted by inhibition
- `probe` : Interval between encoding and access
"""
function measure_spike_time(simtime, distal_rate, inh_rate, probe; inh_comp=nothing, samples=200)
    out_spikes = zeros(samples)
    for x in 1:samples
        model = TN.H_distal_proximal
        spikes = protocol_memory(simtime,distal_rate=distal_rate, 
                inh_rate=inh_rate, inh_comp=inh_comp,
                probe=probe, model=model) 
        v = TN.run_tripod(model.ds, spikes, simtime; model...)
        _probe = 10*(probe+ENCODING_TIME)
        first = findfirst(x->x>0, v[1,_probe:end])
        out_spikes[x] = isnothing(first) ? simtime : (first*0.1)
    end
    return out_spikes #, mean(spikes)
end


"""
Run simulation for memory encoding and access.
"""

function test_memory(inputs)
    @unpack simtime, inh_comp, inh_rate, probes, distal_rate= inputs
    @info inputs
    samples = 300
    spikes = zeros(samples,length(distal_rate), length(inh_rate), length(probes))
    for p in ProgressBar(eachindex(probes))
    Threads.@threads for n in eachindex(distal_rate)
            for m in eachindex(inh_rate)
                _s = measure_spike_time(simtime, distal_rate[n], inh_rate[m], probes[p]; inh_comp=inh_comp, samples=samples)
                spikes[:,n,m, p] = _s
            end
        end
    end
    return spikes
end

k = 1000
BASE_SOMA_RATE = 100.f0
ENCODING_TIME = 50
simtime= 500
probe_time =40
simtime = 500
base_values = (simtime=simtime, 
            inh_comp="", 
            inh_rate=0:0.2k:5k, 
            probes=0:5:150, 
            distal_rate=0:0.2k:5k )

##
base_inputs = Dict(pairs(base_values))
base_inputs[:inh_comp] = ""
base_inputs[:inh_rate] = [0.f0]
no_inh = test_memory((;base_inputs...))

#
dend_inputs = Dict(pairs(base_values))
dend_inputs[:inh_comp] = "dend"
dend_inputs[:probes] = [0,25,50]
dend = test_memory((;dend_inputs...))

#
soma_inputs = Dict(pairs(base_values))
soma_inputs[:inh_comp] = "soma"
soma_inputs[:probes] = [0,25,50]
soma = test_memory((;soma_inputs...))

file =datadir("dendritic_memory","spike_time_data_balanced.jld")
data = @strdict no_inh dend soma base_values
safesave(file, data)