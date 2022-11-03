"""
Measure the of distributing excitation and inhibition on different compartments.
Measure axial currents and membrane potential in the on-path, on-soma, off-path conditions.
The file reproduce figure 6B.
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using StatsPlots

simtime = 10000
νexc = 4.
νs = 0.0:0.4:8

"""
Measure average membrane potential and average axial current while varying the location and the strength of the inhibitory inputs
"""
function test_inhibition(model, comp, νs)
	νinh = [0., 0., 0.]
	v = zeros(length(νs),3)
	c = zeros(length(νs),3)
	for (n,νi) in enumerate(νs)
		νinh[comp+1] = νi
		inputs = [0., νexc, 0., νinh...] # inpust as firing rate per compartment
		voltage, current = TN.run_tripod(inputs, simtime, in_curr=true; model...)
		v[n, :] = mean(voltage, dims=2)[:,1]
		c[n, 2:3] = mean(current, dims=2)[:,1]
		c[n,1] = mean(-sum(current[1:2,:],dims=1)[1,:])
	end
	return v, c
end

function inhibition_frequency(model)
	vs_volt = []
	vs_current = []
	for i in [2,0,1]
		voltage,current = test_inhibition(model,i, νs)
		push!(vs_volt,voltage)
		push!(vs_current,current)
	end
	return  vs_volt, vs_current
end

function plot_membrane_current(vs_voltage, vs_current)
	labels= ["soma" "dendrite" "dendrite"]
	lines = [:solid :solid :solid ]
	colors = [:black RED GREEN]
	plots_volt=[]
	plots_current=[]
	for (voltage, current) in zip(vs_voltage,vs_current)
		pv = plot(νs, voltage,  linewidth=3, yticks=false,labels=labels, linestyle=lines,color=colors, legend=false)
		pc = plot(νs,current , color=colors,yticks=false,linewidth=3, label="")
		push!(plots_volt,pv)
		push!(plots_current,pc)
	end
	plot!(plots_volt[1], ylabel="Membrane (mV)", yticks=true, legend=true, )
	p =plot(plots_volt..., layout =(1,3), frame=:axes, ylims=(-90,0))
	plot!(plots_current[1], ylabel="Currents (pA)", legend=true,  yticks=(-2000:2000:2000, -2:2:2), ylims=(-3000,3000),xlabel="Inhibition rate (KHz)")
	q = plot(plots_current..., layout=(1,3), frame=:axes)
	return p,q
end

##
model = TN.H_distal_distal
mem, curr = inhibition_frequency(model)
p,q = plot_membrane_current(mem, curr)
z = plot(p,q, layout=(2,1))

savefig(z,plotsdir("dendritic_memory","Fig6B.pdf"))
