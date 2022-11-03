"""
Show the correlation between membrane potential and switch events in the  localization of input. 

First part:
    1. Generate a signal -> oscillating / exp_decay.
	2. Apply the signal to the spikes -> apply_stimuli!
	3. Run the tripod function
	4. Make plots
Second part:
	1. Run same protocl as above, but for 1000 cells.
	2. Measure the average membrane potential.
	4. Make plots
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include("stimuli.jl")
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using HDF5, Random, Measures

##
Random.seed!(45)
simtime = 2_500
signal = oscillating_signal(simtime, 200 )
stimuli = get_stimuli(signal)
inputs = TN.make_spikes(simtime; rate=1200)
apply_stimuli!(stimuli, inputs, soma=false)


default(palette=:blues)
voltage = TN.run_tripod((400,400),inputs,simtime, do_spikes=true,  soma_only=false)
xx = 8_001:18_000
p = plot(voltage[3,xx] , lw=3)
plot!(voltage[2,xx] , lw=3)
plot!(voltage[1,xx], lw=2, c=:black, alpha=0.7)
plot!(legend=false, xticks=(range(1,10_000,5), 0:0.25:1))
plot!(guidefontsize=18, tickfontsize=13, xlabel="Time (ms)", grid=false, yticks=false)
for x in 0:1:4
	plot!([(x)*2000,(1+x)*2000],[-80, -80], lw=5,)
end

voltage = TN.run_tripod((400,150),inputs,simtime, do_spikes=true,  soma_only=false)
q = plot(voltage[3,xx] , lw=3)
plot!(voltage[2,xx] , lw=3)
plot!(voltage[1,xx], lw=2, c=:black, alpha=0.7)
plot!(legend=false, xticks=(range(1,10_000,5), 0:0.25:1))
plot!(guidefontsize=18, tickfontsize=13, xlabel="Time (ms)", grid=false, yticks=false)
for x in 0:1:4
	plot!([(x)*2000,(1+x)*2000],[-80, -80], lw=4,)
end
plot!(p, ylabel="Membrane potential (mV)")
##
"""
Average the somatic potential over 2000 different realizations of the 
same input protocol.
"""
file = datadir("dendritic_memory","membrane_potential_transition.jld")

if isfile(file)
	mean_v = JLD.load(file)["mean_v"]
else
	models = TN.models[[1,3]]
	samples = 2000
	simtime = 2500
	mean_v = zeros(10_000,samples,2);
	stimuli = get_stimuli(signal);
	AdEx = TN.AdExParams(up=0.1, idle=0.0)
	Threads.@threads for x in 1:samples
		inputs = TN.make_spikes(simtime; rate=1200)
		inh_inputs = TN.make_spikes(simtime; rate=1000)
		inputs[4:6,:] = inh_inputs[4:6,:]
		apply_stimuli!(stimuli, inputs)
		for m in eachindex(models)
			voltage = TN.run_tripod(models[m].ds,inputs,simtime, do_spikes=true,  soma_only=false)
			mean_v[:,x,m] .=  voltage[1,12_001:22_000]
		end
	end

	data = @strdict mean_v
	safesave(file, data)
end

average_membrane = mean(mean_v, dims=2)[:,1,:] .-mean(mean_v)
plot!(p, 10*(average_membrane[:,1] ) .-140, c=:black, lw=4)
plot!(q, 7*(average_membrane[:,2] ) .-120, c=:black, lw=4)

z = plot(p,q, margins= 5mm, ylims=(-150,20))
savefig(z, plotsdir("dendritic_memory","Fig9B.pdf"))
z