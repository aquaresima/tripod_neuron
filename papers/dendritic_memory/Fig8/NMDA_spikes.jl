"""
Compare EPSP duration between the synapses with mouse parameters and synapses with human parameters.
After, measure the duration of the dendritic EPSP
This file reproduces Fig 8A and Fig 8B
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

##

function single_spike(tt::Int64, comp::Int64 ;eff=1.f0)
    comp == 2 && tt == EXCSPIKETIME && return eff
    return 0
end


DEPOLARIZED_THRESHOLD = -60
function NMDA_spikes( ;delay=0, model, my_range = 0:15:200)
    s=plot()
    p=plot()
    above_threshold = Vector()
    end_threshold = Vector()
	CList = palette(:roma, length(my_range) );

    for (n,intensity) in enumerate(my_range)
        voltage = TN.run_tripod(model.ds, single_spike; model..., eff=intensity)

        # Measure all the points above the 
	    push!(above_threshold, TN.dt*length(findall(x->x>= DEPOLARIZED_THRESHOLD, voltage[2,EXCSPIKETIME:end])))
	    push!(end_threshold, findlast(x->x>DEPOLARIZED_THRESHOLD, voltage[2,1:end]))

	    s= plot!(s, [voltage[2,:] ],label="", color=CList[n], linewidth=2)
	    p= plot!(p, [voltage[1,:] ], ylims=(-77, -50), color=CList[n], label="", linewidth=2)
	end

	plot!(p,ylabel="Soma (mV)")
	plot!(s,ylabel="Dendrite (mV)")

    s =xlabel!(s, "")
    xlims!(p,(1,2000))
    xlims!(s,(1,2000))
    p = xticks!(p,1001:500:2700, string.(0:50:170))
    end_threshold[1] =1000
    end_threshold[2] =1100
    return s,p, above_threshold
end

## Run

EXCSPIKETIME = 150
model = (ds = (400,400), syn_model = TN.human_synapses, species="H")
s1,p1, aboveNMDA= NMDA_spikes(model=model)

model = (ds = (400,400), syn_model = TN.mouse_synapses, species="H")
s,p, aboveAMPA = NMDA_spikes(model =model);

plot!(p, xlabel="Time after spike (ms)")
p1 = plot(p1, frame=:axes, ylabel="");
s1 = plot(s1, frame=:axes, ylabel="");
plot!(p1, yticks=nothing)
plot!(s1, ticks=nothing, title="Human (NAR = 1.8)");
plot!(s, xticks=nothing ,title="Mouse (NAR = 0.25)");
##
z = plot(s,s1,p,p1, layout=(2,2),frame=:axes )
hm = hcat([[n] for n in 0:10]...)
inset = (1, bbox(0.15, 0.55, 0.6, 0.1, :bottom, :right))
heatmap!(0:20:200, [1], hm,c=:roma, inset=inset,
                    subplot=5, cbar=false, 
                    yaxis=false, 
                    yticks=:none,rotation=-45)

savefig(z,plotsdir("dendritic_memory","Fig8A.pdf"))

"""
Measure the duration of the plateau potential for different dendritic lengths and efficacies
"""
function above_threshold(ls,effs)
    plateau = plot()
	soma_spikes = zeros(length(ls))
    delay=0
    plateau = zeros(length(ls),length(effs))
    for n in eachindex(ls)
		l = ls[n]
        model = (ds = (l,l), syn_model = TN.human_synapses, species="H")
        _,_, plateau_NMDA = NMDA_spikes(delay=delay,  model=model, my_range=effs)
        plateau[ n, :] =  plateau_NMDA
    end

	## make plots
	for l in 1:length(ls)
		for x in 1:length(effs)-1
			if plateau[l,x+1] < plateau[l,x]
				soma_spikes[l] = x
				# plateau[1,l,x:end] .= plateau[1,l,x]
				break
			end
		end
	end
	return plateau, soma_spikes
end

effs = 1:2:200
ls = 100:30:500
plateau, soma_spikes = above_threshold(ls,effs)
##
cList = palette(:blues, length(ls), rev=true)
p = plot(effs, plateau', ylims=(0,120), lw=2, c=[c for c in cList]', label="")
for n in eachindex(soma_spikes)
	if soma_spikes[n]>0
		s = Int(soma_spikes[n])
		plot!(effs[s:end],plateau[n,s:end], c=:white, ylims=(0,120), lw=3, alpha=0.8, label="")
		scatter!([effs[s]],plateau[n,s:s], m=:diamond, ms=8, c= cList[n], msc=cList[n], labels="")
	end
end
plot!(title="NMDA plateau potential", titlefontsize=18,
        ylabel="Plateau duration (ms)", xlabel="No coincident pre-synaptic spikes")
myfont = Plots.text("").font
myfont.rotation=14
plot!(ylims=(0,125),guidefontsize=18,tickfontsize=13, frame=:axes)
inset = (1, bbox(0.6, 0.8, 0.30, 0.07, :bottom, :right))
z = hcat([[n] for n in reverse(1:length(ls))]...)
heatmap!(ls,[1],z,c=:blues, inset=inset, subplot=2, cbar=false, yticks=:none,
                        rotation=-45, yaxis=false, frame=:axes, title="dendritic length (μm)")

savefig(p,plotsdir("dendritic_memory","Fig8B.pdf"))
##