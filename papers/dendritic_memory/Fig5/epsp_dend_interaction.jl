# Copyright (c) 2022 Author Name (AQ)
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


"""
This file measures the effect of interaction between two excitatory conductances and reproduces Fig 4a and Fig 4b
"""

using LsqFit,IterTools, StatsBase

"""
Protocol with two inputs on same or different dendritic compartment
"""
function double_input(tt::Int64, comp::Int; g1::Real, g2::Real, spiketime::Int64, same=false)
	if tt==spiketime
		if same 
			comp==2 && return g1+g2
		else
			comp==2 && return g1
			comp==3 && return g2
		end
	end
	return 0
end



"""
Run the simulation setting the NMDA receptor, the conductance and the model (dendritic lengths)
"""
function gg_simulate(g1::Real,g2::Real; model, same=false, NMDA=false)
	model = Dict(pairs(model))
	if NMDA
		model[:syn_model] =TN.default_synapses
	else
		model[:syn_model] =TN.ampa_only
	end
	model = (; model...)
	spiketime = 1000
	voltage = TN.run_tripod(model.ds,double_input, 500; 
					model..., spiketime=spiketime, same=same,
					g1=g1, g2=g2)
	epsp = TN.get_PSP_peak(voltage, spiketime=spiketime)
    return epsp

end


function gg_synaptic_efficacy(;gs=nothing, func::Function, kwargs...)
	Δ = Vector()
	g = Vector()
	for (a,b) in eachrow(gs)
		push!(Δ, func(a,b;kwargs...))
		push!(g,a*b)
	end
	return Δ, g
end
