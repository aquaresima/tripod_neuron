##
# include("../STDP_params.jl")
using ProgressBars
struct ISTDP
    τy::Float32
    η::Float32
    r0::Float32
    α::Float32
    j⁻::Float32  # minimum weight
    j⁺::Float32  # maximum weight
end

function run_tripod_vogels(
    model::Tuple{Int,Int},
    inputs::Union{Matrix{Bool},Matrix{Int},Vector{Float32}},
    simtime::Int,
)

    pm1 = PassiveMembraneParameters("d1", "H", 4, model[1])
    pm2 = PassiveMembraneParameters("d2", "H", 4, model[2])
    ## Define variables
    v_soma = AdEx.Er
    w_soma = 0.0

    c_soma = zeros(Float32, 2) # soma currents
    v_dend1 = AdEx.Er
    v_dend2 = AdEx.Er

    # Spike events
    idle_soma = 1 # time idle for soma
    soma_th = AdEx.θ  # threshold for spike

    ## structured like:
    # (h,g) x (AMPA, NMDA, GABAa, GABAb)
    syn_soma = zeros(Float32, 4, 2)
    syn_dend1 = zeros(Float32, 4, 2)
    syn_dend2 = zeros(Float32, 4, 2)

    ## Inputs
    exc_soma = 0.0
    exc_dend1 = 0.0
    exc_dend2 = 0.0
    inh_soma = 0.0
    inh_dend1 = 0.0
    inh_dend2 = 0.0


    spiked_exc = false
    ##

    #Set simulation
    total_steps = round(Int, (simtime) / dt)
    iterations = ProgressBar(1:total_steps)
    v = zeros(3, total_steps) # voltage trace
    # e = zeros(total_steps,3) # voltage trace

    ## Duplet traces update before learning rule
    rate = 0.0
    i_trace = zeros(3)
    inh_prespike = falses(3)
    weights = zeros(3, total_steps)
    GABAsynapses = ones(3)
    ##


    # for tt in iterations
    for tt = 1:total_steps
        # Get noise inpute
        if isa(inputs, Vector{Float32})
            exc_soma = 1.0 * rand(Poisson(inputs[1] * dt))
            exc_dend1 = 1.0 * rand(Poisson(inputs[2] * dt))
            exc_dend2 = 1.0 * rand(Poisson(inputs[3] * dt))
            inh_soma = GABAsynapses[1] * rand(Poisson(inputs[4] * dt))
            inh_dend1 = GABAsynapses[2] * rand(Poisson(inputs[5] * dt))
            inh_dend2 = GABAsynapses[3] * rand(Poisson(inputs[6] * dt))
            # @show inh_dend1
        end


        if isa(inputs, Matrix{Bool})
            exc_soma = 1.0 * inputs[1, tt]
            exc_dend1 = 1.0 * inputs[2, tt]
            exc_dend2 = 1.0 * inputs[3, tt]
            inh_soma = 1.0 * inputs[4, tt]
            inh_dend1 = 1.0 * inputs[5, tt]
            inh_dend2 = 1.0 * inputs[6, tt]
        end
        # Get prepared input

        if isa(inputs, Matrix{Int})
            exc_soma = 1.0 * inputs[1, tt]
            exc_dend1 = 1.0 * inputs[2, tt]
            exc_dend2 = 1.0 * inputs[3, tt]
            inh_soma = 1.0 * inputs[4, tt]
            inh_dend1 = 1.0 * inputs[5, tt]
            inh_dend2 = 1.0 * inputs[6, tt]
        end
        # e[tt,:] = [exc_soma, exc_dend1, exc_dend2]
        #
        @views update_synapse_soma!(syn_soma[:, :], inh_soma, exc_soma, Esyn_soma)
        @views update_synapse_dend!(syn_dend1[:, :], inh_dend1, exc_dend1, Esyn_dend)
        @views update_synapse_dend!(syn_dend2[:, :], inh_dend2, exc_dend2, Esyn_dend)


        ## Threshold adaptation decay
        soma_th -= dt * postspike.τA⁻ * (soma_th - AdEx.θ)

        @assert(!isnan(v_soma))

        idle_soma -= 1

        if idle_soma < 0
            ## compute currents
            c_soma[1] = -(v_dend1 - v_soma) * pm1.g_ax
            c_soma[2] = -(v_dend2 - v_soma) * pm2.g_ax
            ## apply currents
            v_dend1 += dt * pm1.C⁻ * c_soma[1]
            v_dend2 += dt * pm2.C⁻ * c_soma[2]
            v_soma -= dt * AdEx.C⁻ * (c_soma[1] + c_soma[2])
            ## update dendritic compartments
            v_dend1 +=
                v_dend1 > 0.0 ? 0.0 : dt * ΔvDend(v_dend1, (@view syn_dend1[:, 2]), pm1)
            v_dend2 +=
                v_dend2 > 0.0 ? 0.0 : dt * ΔvDend(v_dend2, (@view syn_dend2[:, 2]), pm2)
            # update soma
            ## Spike behaviour
            if v_soma >= AdEx.θ + 10
                spiked_exc = true
                ## Increase spike thresshold and set the idle period
                idle_soma = AdEx.up_dt + AdEx.idle_dt
                v_soma = AdEx.AP_membrane
                soma_th += postspike.A
                w_soma += AdEx.b
            else
                v_soma += dt * Δv(v_soma, w_soma, soma_th, (@view syn_soma[:, 2]))
                w_soma += dt * Δw(v_soma, w_soma, soma_th, (@view syn_soma[:, 2]))
            end
        elseif idle_soma > AdEx.idle_dt
            v_soma = AdEx.BAP
            c_soma[1] = -(v_dend1 - v_soma) * pm1.g_ax * BAP_gax
            c_soma[2] = -(v_dend2 - v_soma) * pm2.g_ax * BAP_gax
            ## apply currents
            v_dend1 += dt * pm1.C⁻ * c_soma[1]
            v_dend2 += dt * pm2.C⁻ * c_soma[2]
            v_soma -= dt * AdEx.C⁻ * (c_soma[1] + c_soma[2])
            ## update dendritic compartments
            v_dend1 +=
                v_dend1 > 0.0 ? 0.0 : dt * ΔvDend(v_dend1, (@view syn_dend1[:, 2]), pm1)
            v_dend2 +=
                v_dend2 > 0.0 ? 0.0 : dt * ΔvDend(v_dend2, (@view syn_dend2[:, 2]), pm2)
        else
            v_soma = AdEx.Er
        end


        inh_prespike[1] = inh_soma > 0
        inh_prespike[2] = inh_dend1 > 0
        inh_prespike[3] = inh_dend2 > 0

        # @show inh_prespike
        rate -= dt * rate / iSTDP.τy
        i_trace += inh_prespike - dt * i_trace / iSTDP.τy

        if spiked_exc
            rate += 1.0
        end
        for d = 1:3
            if inh_prespike[d] && GABAsynapses[d] >= 0.0
                GABAsynapses[d] += iSTDP.η * (rate - iSTDP.α)
            end
            if spiked_exc
                GABAsynapses[d] += iSTDP.η * i_trace[d]
            end
        end
        spiked_exc = false

        weights[:, tt] .= GABAsynapses

        v[1, tt] = v_soma
        v[2, tt] = v_dend1
        v[3, tt] = v_dend2
    end
    return v, weights, rate
end
