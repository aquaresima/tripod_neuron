using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
include("epsp_dend_interaction.jl")


"""
Measure the EPSP_peak difference between the AA and the AB condition
"""
function clustered_spread_interaction(a,b ; kwargs...)
	AA = gg_simulate(a,b;same=true, kwargs...)
	AB = gg_simulate(a,b;same=false, kwargs...)
	return AA - AB
end


function clustered_spread(models)
	gs = rand(1.:35,200,2)
	results = zeros(length(models), 2, size(gs)[1]) # model, NMDA/AMPA, AA/AB, gs
	conditions = [true, false]
	Threads.@threads  for n in eachindex(models)
		model= models[n]
		for (nmda,NMDA) in enumerate(conditions)
			Δ, _ = gg_synaptic_efficacy(gs=gs, func=clustered_spread_interaction, NMDA=NMDA, model=model )
			results[n, nmda, :] .= Δ
		end
	end
	return results, gs
end

## Compare AA and AB
models = TN.models[[1,2]]
Δepsp, synapses = clustered_spread(models)
gg_sum = sum(synapses,dims=2)
plots = []
for model in 1:2
	p = Plots.scatter(gg_sum, Δepsp[model,1,:] ,markersize=6, c=RED, msc=RED);
	Plots.scatter!(gg_sum,Δepsp[model,2,:],markersize=6, c=BLU, msc=BLU);
	Plots.plot!([minimum(gg_sum), maximum(gg_sum)],[0,0], ls=:dash, lw=3, c=:black)
	Plots.plot!(xlims=(minimum(gg_sum), maximum(gg_sum)), ylims=(-5,5), legend=false)
	if model ==1
		plot!(title="distal-distal", ylabel=L"\Delta"*"EPSP (mV)")
	end
	if model ==2
		plot!(xlabel="Total "*L" |g_e| (nS)", title="proximal-proximal")
	end
	push!(plots,p)
end
p = plot(plots[[2,1]]..., layout=(1,2), ylims=(-5,4), 
                xlims= (0,85),yticks=([-4,0,3],[-4,0,3]))

savefig(p,plotsdir("dendritic_memory","Fig5A.pdf"))

