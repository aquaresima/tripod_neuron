# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


#===================================================
		Backpropagation spike
===================================================#

# Set AP_membrane for backpropagation
@with_kw struct PostSpike
		# After spike adaptive threshold
		A::Float32 = 10.
		τA::Float32 = 30.
 		Ips::Float32=0. # Introduced in Clopath for vSTDP
		# After spike timescales and membrane
		AP_membrane::Float32=21f0
		BAP::Float32=20f0
		τA⁻::Float32 = 1/τA
		τz⁻::Float32=1.
        BAP_gax::Float32=1.f0
end



#===================================================
                Neuron struct
===================================================#

abstract type NeuronParams end

@with_kw struct AdExParams <: NeuronParams
    #Membrane parameters
    C::Float32=281 # (pF) membrane timescale
    gl::Float32=40   # (nS) gl is the leaking conductance,opposite of Rm
    Rm::Float32=1/gl# (GΩ) total membrane resistance
    τm::Float32=C/gl # (ms) C / gl
    Er::Float32=-70.6 # (mV) resting potential

    # AdEx model
    u_r::Float32=-70.6 # (mV) Reset potential of membrane
    θ::Float32=-50.4 # (mv) Rheobase threshold
    ΔT::Float32=2 # (mV) Threshold sharpness

    # Adaptation parameters
    τw::Float32=144 #ms adaptation current relaxing time
    a::Float32=4 #nS adaptation current to membrane
    b::Float32=80.5 #pA adaptation current increase due to spike

	up::Float32=1 #ms
	idle::Float32=2 #ms

    # Inverse value for simulation speedup
    C⁻::Float32=1/C  # (pF) inverse membrane timescale
    τw⁻::Float32=1/τw #ms inverse adaptation current relaxing time
    τm⁻::Float32=1/τm #ms inverse adaptation current relaxing time
    ΔT⁻::Float32=1/ΔT # (mV) inverse Threshold sharpness
end
