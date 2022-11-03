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


##
k = 1000
BASE_SOMA_RATE = 100.f0
ENCODING_TIME = 50
simtime= 500
probe_time =20
simtime = 500

base_values = (simtime=simtime, 
            inh_comp="", 
            inh_rate=0:0.2k:5k, 
            probes=0:5:150, 
            distal_rate=0:0.2k:5k )

model = TN.H_distal_proximal
spikes = protocol_memory(simtime, distal_rate=1k, inh_rate=0k, 
                        inh_comp="dend", probe=probe_time, model=model) 

##
v,w = TN.run_tripod(model.ds, spikes, simtime, adapt=true;)
default(lw=1)
p = plot(v[1,:], label="soma", c=:black)
plot!(v[2,:], label="Distal dendrite")
plot!(v[3,:], label="Proximal dendrite")
plot!(xticks=(range((ENCODING_TIME+probe_time)*10,5000,5), range(0,400,5)), ylims=(-90,90), xlabel="Time (ms)", ylabel="Membrane potential (mV)")
annotate!(-1235, 140, Plots.text("C", 18, :black, :bold))
plot!(ylims=(-90,150), legendfontsize=10)
plot!(yticks=(-80:20:20, -80:20:20))
q = twinx(p)
plot!(q,w, xticks=:none, ylabel="Adapt. current (pA)", label="")
vline!([500], label="End encoding phase", c=:black, ls=:dash, lw=3)
vline!([500+probe_time*10], label="Retrieval phase", c=:red, ls=:dash, lw=3)
plot!(p,right_margin=20Plots.mm, legend_background_color=:white, legend=:topleft)
plot!(size=(600,600), legend_bg=:white)
savefig(p,plotsdir("dendritic_memory","Fig8_appendix.pdf"))
p