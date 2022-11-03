include("equations.jl")

run_tripod(inputs,simtime; kwargs...) = run_tripod((-1,-1),inputs,simtime; kwargs...)

function run_tripod(model::Union{Tuple{Int, Int}, Vector{Int}}=(-1,-1), 
	inputs::Union{Matrix{Bool}, Matrix{Int}, Matrix{Float32}, Vector{<:Real}, Function}=null_input, 
	simtime::Int=1000;
	species::String="H",
	ds::Union{Tuple{Int, Int}, Vector{Int}} = (-1, -1),
	syn_model::TripodSynapses=default_synapses,
	AdEx::AdExParams=AdExTripod,
	postspike=PostSpike(),
	synapses=[1.,1.,1,1.,1.,1.], 
	ext_currents=zeros(3),
	adapt=false,
	soma_only=false,
	do_spikes=true,
	in_curr = false,
	rec_syn=false,
	kwargs...
	)

	c_soma= zeros(Float32,2) # soma currents
	@unpack Esyn_dend, Esyn_soma = syn_model
	idle_up_dt= AdEx.up/dt + AdEx.idle/dt
	idle_dt=  AdEx.idle/dt

	## Compatibility with old syntax
	if sum(model )<0
		model = ds
	end
	if sum(model) == 0 
		soma_only || @warn "Soma only not set correctly"
		soma_only=true
	end
	@assert sum(model) > -1
	if species=="human"
		species="H"
	end
	if species=="mouse"
		species="M"
	end
	@debug "Model: ", model, "Simtime: ", simtime, "Species: ", species,  "Adapt: ", adapt, "Soma only: ", soma_only, "Do spikes: ", do_spikes, "In curr: ", in_curr
	##

	pm1   = PassiveMembraneParameters("d1", species, 4, model[1])
	pm2   = PassiveMembraneParameters("d2", species, 4, model[2])

	@inline ΔvSoma(v,w,axial, θ, g)= ΔvSoma2(v,w,axial, θ, g,AdEx, Esyn_soma)
	@inline ΔvSoma_nospike(v,w,axial, g)= Δv_soma_nospike(v,w,axial, g, AdEx, Esyn_soma)
	@inline ΔwSoma(v,w)= Δw_soma(v,w,AdEx)
	@inline ΔvDend(v,axial, g,  pm) = ΔvDend2(v,axial, g, pm, Esyn_dend)

	## Define variables
	v_soma::Float32    = AdEx.Er
	w_soma::Float32 =    0.
	v_dend1::Float32 =    AdEx.Er
	v_dend2::Float32 =    AdEx.Er
	c_soma= zeros(Float32,2) # soma currents

	# Spike events
	last_event = 0 # time idle for soma
	soma_th = AdEx.θ  # threshold for spike

	## structured like:
	# (h,g) x (AMPA, NMDA, GABAa, GABAb)
	syn_soma =  zeros(Float32,4,2)
	syn_dend1 = zeros(Float32,4,2)
	syn_dend2 = zeros(Float32,4,2)

	## Inputs
	exc_soma::Float32 =  0.f0
	exc_dend1::Float32 = 0.f0
	exc_dend2::Float32 = 0.f0
	inh_soma::Float32 =  0.f0
	inh_dend1::Float32 = 0.f0
	inh_dend2::Float32 = 0.f0


	#Set simulation
	total_steps = round(Int,(simtime)/dt)
	v = zeros(Float32,4,total_steps) # voltage trace
	in_curr && (c = zeros(Float32,2,total_steps) )# voltage trace
	rec_syn && (syn = zeros(Float32,3,total_steps, 4)) # voltage trace

	# for tt in iterations
	for tt in 1:total_steps

		# Get external inputs
	    v_soma  += dt*AdEx.C⁻*ext_currents[1]
		v_dend1 += dt*pm1.C⁻* ext_currents[2]
		v_dend2 += dt*pm2.C⁻* ext_currents[3]

		if isa(inputs, Vector)
			exc_soma  = synapses[1] *rand(Poisson(inputs[1]*dt))
			exc_dend1 = synapses[2] *rand(Poisson(inputs[2]*dt))
			exc_dend2 = synapses[3] *rand(Poisson(inputs[3]*dt))
			inh_soma  = synapses[4] *rand(Poisson(inputs[4]*dt))
			inh_dend1 = synapses[5] *rand(Poisson(inputs[5]*dt))
			inh_dend2 = synapses[6] *rand(Poisson(inputs[6]*dt))
		elseif isa(inputs, Matrix{Bool})
			exc_soma  =synapses[1] * inputs[1,tt]
			exc_dend1 =synapses[2] * inputs[2,tt]
			exc_dend2 =synapses[3] * inputs[3,tt]
			inh_soma  =synapses[4] * inputs[4,tt]
			inh_dend1 =synapses[5] * inputs[5,tt]
			inh_dend2 =synapses[6] * inputs[6,tt]
		elseif isa(inputs, Matrix{Float32}) || isa(inputs, Matrix{Int})
			exc_soma  =synapses[1] * inputs[1,tt]
			exc_dend1 =synapses[2] * inputs[2,tt]
			exc_dend2 =synapses[3] * inputs[3,tt]
			inh_soma  =synapses[4] * inputs[4,tt]
			inh_dend1 =synapses[5] * inputs[5,tt]
			inh_dend2 =synapses[6] * inputs[6,tt]
		elseif isa(inputs, Function)
			exc_soma = inputs(tt,1; kwargs...)
			exc_dend1= inputs(tt,2; kwargs...)
			exc_dend2= inputs(tt,3; kwargs...)
			inh_soma = inputs(tt,4; kwargs...)
			inh_dend1= inputs(tt,5; kwargs...)
			inh_dend2= inputs(tt,6; kwargs...)
		end

		if soma_only
			# @views update_synapse_soma!(syn_soma[:,:],   inh_soma, exc_soma, Esyn_soma)
			## Treat each of the somatic synapses as separate, 
			## which is equivalent to dendrites of zero length. 
			@views update_synapse_soma!(syn_soma[:,:],   inh_soma, exc_soma, Esyn_soma)
			@views update_synapse_soma!(syn_dend1[:,:],  inh_dend1, exc_dend1, Esyn_soma)
			@views update_synapse_soma!(syn_dend2[:,:],  inh_dend2, exc_dend2, Esyn_soma)
			syn_soma = .+ syn_dend1 .+syn_dend2
		else 
			@views update_synapse_soma!(syn_soma[:,:],   inh_soma, exc_soma, Esyn_soma)
			@views update_synapse_dend!( syn_dend1[:,:], inh_dend1, exc_dend1, Esyn_dend)
			@views update_synapse_dend!( syn_dend2[:,:], inh_dend2, exc_dend2, Esyn_dend)
		end
		if rec_syn
			Rs = fieldnames(TN.Synapse)
			syn[1,tt,:] = [g* getfield(Esyn_soma,R).gsyn for (g,R) in zip(syn_soma[:,2], Rs)]
			syn[2,tt,:] = [g* getfield(Esyn_dend,R).gsyn for (g,R) in zip(syn_dend1[:,2], Rs)]
			syn[3,tt,:] = syn_dend2[:,2]
		end


		@assert(!isnan(v_soma))
		## Diminish last_event counter
        last_event -= 1
		# Precompute currents before the voltage update
		if soma_only
			fill!(c_soma,0.f0)
		else
			c_soma[1] = - (v_dend1 - v_soma)* pm1.g_ax
			c_soma[2] = - (v_dend2 - v_soma)* pm2.g_ax
		end
		# ordinary neuronal integration: no spike, no reset 
		if last_event < 0
			# check if super-threshold
			if !do_spikes
				v_dend1 += dt*ΔvDend(v_dend1, -c_soma[1],(@view syn_dend1[:,2]), pm1)
				v_dend2 += dt*ΔvDend(v_dend2, -c_soma[2],(@view syn_dend2[:,2]), pm2)
				v_soma +=  dt*ΔvSoma_nospike(v_soma, w_soma, sum(c_soma), (@view syn_soma[:,2]))
			else
				# spike if above threshold
				if v_soma >= AdEx.θ
					spiked_exc = true
					last_event = idle_up_dt 
					v_soma  = postspike.AP_membrane
					w_soma  += AdEx.b
				# otherwise integrate
				else
					v_dend1 += dt*ΔvDend(v_dend1, -c_soma[1],(@view syn_dend1[:,2]), pm1)
					v_dend2 += dt*ΔvDend(v_dend2, -c_soma[2],(@view syn_dend2[:,2]), pm2)
					w_soma += dt*ΔwSoma(v_soma, w_soma)
					v_soma += dt*ΔvSoma(v_soma, w_soma, sum(c_soma), soma_th, (@view syn_soma[:,2]))
				end
			end
		## Action Potential driven backpropagation
		elseif last_event > idle_dt
	        v_soma = postspike.BAP
			v_dend1 += dt*ΔvDend(v_dend1, -c_soma[1],(@view syn_dend1[:,2]), pm1)*postspike.BAP_gax
			v_dend2 += dt*ΔvDend(v_dend2, -c_soma[2],(@view syn_dend2[:,2]), pm2)*postspike.BAP_gax

		## reset
		else
			v_dend1 += dt*ΔvDend(v_dend1, -c_soma[1],(@view syn_dend1[:,2]), pm1)
			v_dend2 += dt*ΔvDend(v_dend2, -c_soma[2],(@view syn_dend2[:,2]), pm2)
	        v_soma = AdEx.u_r
		end
		if soma_only
			v_dend1  = v_soma
			v_dend2  = v_soma
		end

		## Recordings
		if in_curr
			c[:,tt] = c_soma
			fill!(c_soma,0.f0)
		end


		v[1,tt] = v_soma
		v[2,tt] = v_dend1
		v[3,tt] = v_dend2
		v[4,tt] = w_soma
	end
	if rec_syn
		return syn
	end
	if adapt
		return v[1:3,:], v[4,:]
	end
	if in_curr
		return v[1:3,:], c
	end
	return v[1:3,:]
end
