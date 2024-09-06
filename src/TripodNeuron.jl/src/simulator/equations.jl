@inline function exp32(x::Float32)::Float32
    x = ifelse(x < -10.0f0, -32.0f0, x)
    x = 1.0f0 + x / 32.0f0
    x *= x
    x *= x
    x *= x
    x *= x
    x *= x
    return x
end

@inline function exp256(x::Float32)::Float32
    x = ifelse(x < -10.0f0, -256.0f0, x)
    x = 1.0f0 + x / 256.0f0
    x *= x
    x *= x
    x *= x
    x *= x
    x *= x
    x *= x
    x *= x
    x *= x
    return x
end


@inline function NMDA_nonlinear(NMDA::ReceptorVoltage, v::Float32)::Float32
    Mg_mM = 1.0
    return (1 + (Mg_mM / NMDA.b) * exp32(NMDA.k * (v)))^-1 ##NMDA
end

@inline function syn_current(v::Float32, g::AbstractVector{Float32}, syn::Synapse)::Float32
    return (syn.AMPA.gsyn * g[1] + syn.NMDA.gsyn * g[2] * NMDA_nonlinear(syn.NMDA, v)) *
           (v) +
           syn.GABAa.gsyn * (v - syn.GABAa.E_rev) * g[3] +
           syn.GABAb.gsyn * (v - syn.GABAb.E_rev) * g[4]
end

##  Integration methods
#Heun integration
# y' = f(y)
# γ = yₜ + dt f(y)
# y ₜ₊₁ = yₜ + dt/2 *(f(γ) + f(yₜ)

# # Runge Kutta4 integration
# @inline function ΔvDend4(v::Float32,axial::Float32, g::AbstractVector{Float32},  pm::PassiveMembraneParameters)::Float32
# 	k1 = _ΔvDend(v, axial, g, pm)
# 	k2 = _ΔvDend(v + k1 * dt/2,axial, g, pm)
# 	k3 = _ΔvDend(v + k2 * dt/2,axial, g, pm)
# 	k4 = _ΔvDend(v + k3 * dt  ,axial, g, pm )
#     return (1/6) * (k1 + 2*k2 + 2*k3 + k4)
# end

# Runge Kutta -2 integration
function ΔvDend2(
    v::Float32,
    axial::Float32,
    g::AbstractVector{Float32},
    pm::PassiveMembraneParameters,
    syn::Synapse,
)::Float32
    k1 = _Δv_dend(v, axial, g, pm, syn)
    k2 = _Δv_dend(v + k1 * dt, axial, g, pm, syn)
    return (1 / 2) * (k1 + k2)
end

function ΔvSoma2(v, w, axial, θ, g, Neuron, syn)::Float32
    k1 = _Δv_soma(v, w, axial, θ, g, Neuron, syn)
    k2 = _Δv_soma(v + k1 * dt, w, axial, θ, g, Neuron, syn)
    return (1 / 2) * (k1 + k2)
end

# Euler Integration
@inline function _Δv_dend(
    v::Float32,
    axial::Float32,
    g::AbstractVector{Float32},
    pm::PassiveMembraneParameters,
    syn::Synapse,
)::Float32
    i = syn_current(v, g, syn)
    return pm.C⁻ * (-(v - pm.Er) / pm.Rm - min(abs(i), 1500) * sign(i) - axial)
end

@inline function _Δv_soma(
    v::Float32,
    w::Float32,
    axial::Float32,
    θ::Float32,
    g::AbstractVector{Float32},
    Neuron::NeuronParams,
    syn::Synapse,
)::Float32
    @unpack gl, C⁻, Er, ΔT, ΔT⁻ = Neuron
    return C⁻ *
           (gl * (-v + Er + ΔT * exp32(ΔT⁻ * (v - θ))) - w - syn_current(v, g, syn) - axial) ## external currents
end

@inline function Δw_soma(v::Float32, w::Float32, Neuron::NeuronParams)::Float32
    @unpack Er, a, τw⁻ = Neuron
    return τw⁻ * (a * (v - Er) - w)
end

@inline function Δv_soma_nospike(
    v::Float32,
    w::Float32,
    axial::Float32,
    g::AbstractVector{Float32},
    Neuron::NeuronParams,
    syn::Synapse,
)::Float32
    @unpack gl, C⁻, Er, ΔT, ΔT⁻ = Neuron
    return C⁻ * (gl * (-v + Er) - w - syn_current(v, g, syn) - axial) ## external currents
end

##

function double_exp(syn_arr::AbstractVector{Float32}, syn)
    syn_arr[2] = exp32(-dt * syn.τd⁻) * (syn_arr[2] + dt * syn_arr[1])#),13500)
    syn_arr[1] = exp32(-dt * syn.τr⁻) * (syn_arr[1])
end



function update_synapse_soma!(
    syn_arr::AbstractMatrix{Float32},
    inh_::Float32,
    exc_::Float32,
    syn::Synapse,
)
    syn_arr[1, 1] += exc_ * syn.AMPA.α
    syn_arr[3, 1] += inh_ * syn.GABAa.α
    @views double_exp(syn_arr[1, :], syn.AMPA)
    @views double_exp(syn_arr[3, :], syn.GABAa)
    # @views double_exp_g(syn_arr[2,:], syn.NMDA)
    # @views double_exp_g(syn_arr[4,:], syn.GABAb)
    return nothing
end

function update_synapse_dend!(
    syn_arr::AbstractMatrix{Float32},
    inh_::Float32,
    exc_::Float32,
    syn::Synapse,
)
    syn_arr[1, 1] += exc_ * syn.AMPA.α
    syn_arr[2, 1] += exc_ * syn.NMDA.α
    syn_arr[3, 1] += inh_ * syn.GABAa.α
    syn_arr[4, 1] += inh_ * syn.GABAb.α
    @views double_exp(syn_arr[1, :], syn.AMPA)
    @views double_exp(syn_arr[2, :], syn.NMDA)
    @views double_exp(syn_arr[3, :], syn.GABAa)
    @views double_exp(syn_arr[4, :], syn.GABAb)
    return nothing
end
