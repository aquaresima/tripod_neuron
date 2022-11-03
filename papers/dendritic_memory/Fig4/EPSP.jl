"""
In this file we analyze the EPSP in four conditions:
	HUMAN physiology	HUMAN synapses
	MOUSE physiology    MOUSE synapses
	for all the dendritic lenghts

This file reproduces Fig 2a and Fig 2b
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));


EXCSPIKETIME = 100
function protocol(tt::Int64, comp::Int64 ;eff=2.f0)
	if comp == 2 && tt == EXCSPIKETIME
        return eff
    else
        return 0
    end
end


# Model variations
function get_model(nar, l, species)
	ampa_g0 = 0.73
	if species == "human"
		## Use human synapses
		MyNMDA = TN.ReceptorVoltage(E_rev=0.0, τr = 8, τd=35., g0 = ampa_g0*nar, k= -0.077f0)
		s = "H"
	else
		## Use mouse synapses
		MyNMDA = TN.ReceptorVoltage(E_rev=0.0, τr = 0.99, τd=100., g0 = ampa_g0*nar, k=-0.062f0)
		s ="M"
	end
	MyGLU = TN.Glutamatergic(TN.DuarteGluDend.AMPA, MyNMDA)
	syn_model = TN.TripodSynapses(
		Esyn_soma = TN.Synapse(TN.DuarteGluSoma, TN.MilesGabaSoma),
		Esyn_dend = TN.Synapse(MyGLU, TN.MilesGabaDend)
	)
	return (ds = (l,l), syn_model = syn_model, species="H")
end


"""
Get EPSP for all dendritic lengths
"""
function epsp(nar, eff, l, species)
	model = get_model(nar,l,species)
	voltage = TN.run_tripod((-1,-1), protocol, 200; model..., eff=eff)
	return TN.get_PSP_peak(voltage;spiketime=EXCSPIKETIME,rest=TN.AdEx.Er)
end
# EPSPs = epsp(1.,60, 100, "human", exc)

"""
Return index of maxima in Vector
"""
function get_maxima(data)
    arg_maxima = []
    for x in 2:length(data)-1
        (data[x] > data[x-1]) && (data[x]>data[x+1]) && (push!(arg_maxima,x))
    end
    return arg_maxima
end

"""
Get non-linear upswing position and size
"""
function detect_non_linearity(x,y)
	data = y
	max_value = 0.
	soma_spikes = (x=-10, y=0.)
	nmda_spikes = (x=-10, y=0.)
	second_derivative = diff(diff(data))
	second_derivative[second_derivative .< 0.2] .=0
	arg_max = get_maxima(second_derivative)
	if length(arg_max) == 0
		return soma_spikes, nmda_spikes
	elseif length(arg_max)<2
		arg = arg_max[1]
		epsp = data[arg+2] - data[arg]
		if epsp >10
			soma_spikes = (x=x[arg], y=epsp)
		else
			nmda_spikes = (x=x[arg], y=epsp)
		end
	else
		arg = arg_max[1]
		epsp = data[arg+2] - data[arg]
		nmda_spikes = (x=x[arg], y=epsp)
		arg = arg_max[2]
		epsp = data[arg+2] - data[arg]
		soma_spikes = (x=x[arg], y=epsp)
	end
	return nmda_spikes, soma_spikes
end

# Plot for all length
ls = 100:10:500
nars = sort([collect(0.:0.2:1.6)...,0.218,1.8])
effs = 0.:4:200
EXCSPIKETIME = 100
CList = reverse(range(BLU, stop=RED,length=length(nars)))
#
EPSP_human = zeros(length(ls), length(nars))
EPSP_mouse = zeros(length(ls), length(nars))
Threads.@threads for (i) in eachindex(ls)
	l = ls[i]
	for (j,nar) in enumerate(nars)
		EPSP_human[i,j] = epsp(nar, 60.,l, "human")
		EPSP_mouse[i,j] = epsp(nar, 60, l, "mouse")
	end
end



default(lw=3)
a = plot(ls, EPSP_human, c=CList', legend=false, title="human")
plot!(ls, EPSP_human[:,end], legend=false, title="human", c=:black)
b = plot(ls, EPSP_mouse, c=CList', legend=false, title="mouse")
plot!(ls, EPSP_mouse[:,3], c=:black, legend=false, title="mouse", yticks=:none,
yaxis=false)
p = plot(a,b, ylims=(2,18))
savefig(p,plotsdir("dendritic_memory","Fig4C.pdf"))
p
##

effs = 0.:2:200
EPSP_human = zeros(length(effs), length(nars))
EPSP_mouse = zeros(length(effs), length(nars))
for (i,eff) in enumerate(effs)
	for (j,nar) in enumerate(nars)
		EPSP_human[i,j] = epsp(nar, eff, 300, "human")
		EPSP_mouse[i,j] = epsp(nar, eff, 300, "mouse")
	end
end
# p = plot(
a = plot(effs, EPSP_human, c=CList', legend=false, title="human")
plot!(effs, EPSP_human[:,end], c=:black, legend=false, title="human")
b=plot(effs, EPSP_mouse, c=CList', legend=false, title="mouse")
plot!(effs, EPSP_mouse[:,3], c=:black, legend=false, title="mouse",yaxis=false,yticks=:none)
savefig(p,plotsdir("dendritic_memory","Fig4D.pdf"))
p = plot(a,b)
plot!(ylims=(4,30))

##
effs = 0.:5:200
EPSP_human = zeros(length(effs), length(ls), length(nars))
EPSP_mouse = zeros(length(effs), length(ls), length(nars))
Threads.@threads for i in eachindex(effs)
	eff = effs[i]
	for (j,l) in enumerate(ls)
		for (k,nar) in enumerate(nars)
			EPSP_human[i,j,k] = epsp(nar, eff, l, "human")
			EPSP_mouse[i,j,k] = epsp(nar, eff, l, "mouse")
		end
	end
end

##
p,q = plot(), plot()
for (j,l) in enumerate(ls)
	for (k,nar) in enumerate(nars)
		c = CList[k]
		(k == length(nars)) && (c = RGBA(0,0,0))
		nmda, soma = detect_non_linearity(effs,EPSP_human[:,j,k])
		if nmda.x > 0
			scatter!(p, [nmda.x], [l], markershape=:circle, c=c, msc=c, label="", ms=5.)
		end
		if soma.x > 0
			scatter!(p, [soma.x], [l], markershape=:diamond, c=c, msc=c, label="", ms=7.)
		end
		c = CList[k]
		(k == 4) && (c = RGBA(0,0,0))
		nmda, soma = detect_non_linearity(effs,EPSP_mouse[:,j,k])
		if nmda.x > 0
			scatter!(q, [nmda.x], [l], markershape=:circle, c=c, msc=c, label="", ms=5.)
		end
		if soma.x > 0
			scatter!(q, [soma.x], [l], markershape=:diamond, c=c, msc=c, label="", ms=7.)
		end
	end
end
plot!(q, yticks=:none, yaxis=false)
p = plot!(p,q,legend=false, ylims= (100,520), xlims=(0,220))
savefig(p,plotsdir("dendritic_memory","Fig4B.pdf"))
