##
using LinearAlgebra
using StatsBase
using ProgressBars
using Main


include(joinpath(@__DIR__, "STDP_params.jl"))
ciao = "ciao"
module TripodNeuron
using Logging, Printf, Distributions
include("../TripodNeuron/base.jl")
Base.ciao
@show ciao
# AdEx = get_bonoclopath_params()
# global clopath_vstdp_neuron
# AdEx= clopath_vstdp_neuron()
postspike_vSTDP = PostSpike(-30.4, round(Int, 1 / dt), 1 / 50, 400.0, 1 / 40.0)
end

TN = TripodNeuron

triplet = pfister_visualcortex()
vstdp = clopath_vstdp_visualcortex()

##

function simulate_clopath(
    s::TN.Soma,
    simTime::Int64,
    NMDAsynapses::Vector{Float32},
    AMPAsynapses::Vector{Float32},
    GABAsynapses::Vector{Float32},
    spiketimes::Vector{Int64},
    firing_neurons::Vector{Vector{Int}};
    post_firingtimes::Vector{Int},
    inh_τ = 100.0,
    max_efficacy = 5.0,
    normalize = false,
    current::Float32 = 0.0,
    kwargs...,
)

    ## Initialization and loop variables
    total_steps = round(Int, simTime / TN.dt)
    last_spike = -1
    spiked = (s.v > TN.AdEx.θ)

    ## Learning
    u⁻ = 0.0
    u⁺ = 0.0
    xstdp = zeros(size(AMPAsynapses))
    initial_synapses = sum(AMPAsynapses, dims = 1)[1, :]
    norm = initial_synapses
    delay_v = zeros(round(Int, 1 / TN.dt))

    ## Recording

    voltage = zeros(total_steps)
    trace_xstdp = zeros(size(AMPAsynapses)..., total_steps)
    weights = zeros(size(AMPAsynapses)..., total_steps)
    θ = zeros(total_steps)
    firing_rate = zeros(total_steps)
    # rec_triplet =zeros(4,total_steps)

    ## Input spikes
    reverse!(spiketimes)
    reverse!(firing_neurons)
    next_ft = isempty(spiketimes) ? -1 : pop!(spiketimes)
    next_neurons = next_ft > 1 ? pop!(firing_neurons) : 1

    ## Post firing
    reverse!(post_firingtimes)
    next_post_ft = isempty(post_firingtimes) ? -1 : pop!(post_firingtimes)

    ## Preinitialize
    exc_prespikes = falses(size(AMPAsynapses))
    inh_prespikes = falses(size(GABAsynapses))
    inhibitory = Vector{Int}()
    excitatory = Vector{Int}()

    _firing_rate = 0.0
    I_ps = 0.0
    @fastmath @inbounds for tt in ProgressBar(1:total_steps)
        ## Neuron update
        # if tt - last_spike > 10/TN.dt
        if TN.update_AdEx!(s, spiked)
            spiked = true
            last_spike = tt
            s.θ = TN.postspike_vSTDP.A
            _firing_rate += 1.0
            I_ps = TN.postspike_vSTDP.Ips
        else
            I_ps += -TN.dt * I_ps * TN.postspike_vSTDP.τz⁻
            s.θ += -TN.dt * (s.θ - TN.AdEx.θ) * TN.postspike_vSTDP.τA⁻
            s.v += TN.dt * TN.AdEx.C⁻ * I_ps
            spiked = false
        end
        # else
        # 	s.v = TN.AdEx.AP_membrane-1
        # end

        _firing_rate -= TN.dt * _firing_rate / inh_τ


        ## Force post_firing
        while next_post_ft == tt
            s.v = TN.AdEx.AP_membrane
            if length(post_firingtimes) > 0
                next_post_ft = pop!(post_firingtimes)
            else
                next_post_ft = -1
            end
        end

        ## Input
        while next_ft == tt
            excitatory = next_neurons[next_neurons.>0]
            inhibitory = abs.(next_neurons[next_neurons.<0])
            exc_prespikes[excitatory] .= true
            inh_prespikes[inhibitory] .= true

            for neuron in excitatory
                ## Set post-synaptic spike
                TN.exc_spike_plastic!(
                    s,
                    AMPA_eff = AMPAsynapses[neuron],
                    NMDA_eff = NMDAsynapses[neuron],
                )

            end
            for neuron in inhibitory
                TN.inh_spike!(s, eff = GABAsynapses[neuron])
            end
            ## Pull the queue
            if length(spiketimes) > 0
                next_ft = pop!(spiketimes)
                next_neurons = pop!(firing_neurons)
            else
                next_ft = -1
            end
        end
        s.v += current * TN.AdEx.C⁻ * TN.dt

        #===========================================
        						Learning Clopath
        		===========================================#

        if normalize
            if (tt % 200 == 0)
                norm = initial_synapses / sum(AMPAsynapses)
                AMPAsynapses .*= norm
            end
        end

        xstdp .+= exc_prespikes[:]
        xstdp .+= -TN.dt .* (xstdp) .* vstdp.τx⁻

        delay_v[1:end-1] = delay_v[2:end]
        delay_v[end] = s.v

        ## Duplet traces update before learning rule

        # synaptic plasticity depends on 'pre' synaptic neurons
        #LTD
        for syn = 1:length(AMPAsynapses)
            u⁻ += TN.dt * (delay_v[1] - u⁻) * vstdp.τu⁻ #u
            u⁺ += TN.dt * (delay_v[1] - u⁺) * vstdp.τv⁻ #v
            #LTD
            if exc_prespikes[syn]
                if u⁻ > vstdp.θ⁻
                    AMPAsynapses[syn] -= TN.dt * vstdp.a⁻ * (u⁻ - vstdp.θ⁻)
                end
            end
            #LTP
            if (u⁺ > vstdp.θ⁻) && (s.v > vstdp.θ⁺)
                AMPAsynapses[syn] +=
                    TN.dt .* xstdp[syn] .* vstdp.a⁺ * (s.v - vstdp.θ⁺) * (u⁺ - vstdp.θ⁻)
            end
        end
        trace_xstdp[:, tt] .= xstdp
        AMPAsynapses[findall(AMPAsynapses .> max_efficacy)] .= max_efficacy
        AMPAsynapses[findall(AMPAsynapses .< 0.0)] .= 0

        #===========================================
        						Close learning
        		===========================================#



        ## close simulation step
        exc_prespikes .= false
        inh_prespikes .= false

        ## Recordings
        weights[:, tt] .= AMPAsynapses
        voltage[tt] = s.v
        θ[tt] = s.θ
        firing_rate[tt] = _firing_rate

    end
    ## Current and Voltage plot
    return weights, voltage, θ, firing_rate
    # return voltage, rec_synapses, rec_triplet, weights
end

##

function simulate_pfister(
    s::TN.Soma,
    simTime::Int64,
    NMDAsynapses::Vector{Float32},
    AMPAsynapses::Vector{Float32},
    GABAsynapses::Vector{Float32},
    spiketimes::Vector{Int64},
    firing_neurons::Vector{Vector{Int}};
    post_firingtimes::Vector{Int},
    inh_τ = 100.0,
    max_efficacy = 5.0,
    normalize = false,
    current::Float32 = 0.0,
    kwargs...,
)

    ## Initialization and loop variables
    total_steps = round(Int, simTime / TN.dt)
    last_spike = -1
    spiked = (s.v > TN.AdEx.θ)

    ## Learning
    r1stdp = zeros(Float32, size(AMPAsynapses))
    r2stdp = zeros(Float32, size(AMPAsynapses))
    o1stdp = 0.0
    o2stdp = 0.0
    initial_synapses = sum(AMPAsynapses, dims = 1)[1, :]
    norm = initial_synapses

    ## Recording
    voltage = zeros(total_steps)
    firing_rate = zeros(total_steps)
    weights = zeros(size(AMPAsynapses)..., total_steps)
    θ = zeros(total_steps)
    # rec_triplet =zeros(4,total_steps)

    ## Weights AMPA


    ## Input spikes
    reverse!(spiketimes)
    reverse!(firing_neurons)
    next_ft = isempty(spiketimes) ? -1 : pop!(spiketimes)
    next_neurons = next_ft > 1 ? pop!(firing_neurons) : 1

    ## Post firing
    reverse!(post_firingtimes)
    next_post_ft = isempty(post_firingtimes) ? -1 : pop!(post_firingtimes)

    ## Preinitialize
    exc_prespikes = falses(size(AMPAsynapses))
    inh_prespikes = falses(size(GABAsynapses))
    inhibitory = Vector{Int}()
    excitatory = Vector{Int}()

    _firing_rate = 0.0
    I_ps = 0.0
    @fastmath @inbounds for tt in ProgressBar(1:total_steps)
        # for tt in ProgressBar(1:total_steps)

        if TN.update_AdEx!(s, spiked)
            spiked = true
            last_spike = tt
            s.θ = TN.postspike_vSTDP.A
            _firing_rate += 1.0
            # I_ps = TN.postspike_vSTDP.Ips
        else
            # I_ps += - TN.dt * I_ps*TN.postspike_vSTDP.τz⁻
            s.θ += -TN.dt * (s.θ - TN.AdEx.θ) * TN.postspike_vSTDP.τA⁻
            # s.v  += TN.dt * TN.AdEx.C⁻ * I_ps	*102
            spiked = false
        end

        _firing_rate -= TN.dt * _firing_rate / inh_τ



        ## Force post_firing
        while next_post_ft == tt
            s.v = TN.AdEx.AP_membrane
            if length(post_firingtimes) > 0
                next_post_ft = pop!(post_firingtimes)
            else
                next_post_ft = -1
            end
        end

        ## Input
        while next_ft == tt
            excitatory = next_neurons[next_neurons.>0]
            inhibitory = abs.(next_neurons[next_neurons.<0])
            exc_prespikes[excitatory] .= true
            inh_prespikes[inhibitory] .= true

            for neuron in excitatory
                ## Set post-synaptic spike
                TN.exc_spike_plastic!(
                    s,
                    AMPA_eff = AMPAsynapses[neuron],
                    NMDA_eff = NMDAsynapses[neuron],
                )

            end
            for neuron in inhibitory
                TN.inh_spike!(s, eff = GABAsynapses[neuron])
            end
            ## Pull the queue
            if length(spiketimes) > 0
                next_ft = pop!(spiketimes)
                next_neurons = pop!(firing_neurons)
            else
                next_ft = -1
            end
        end

        s.v += current * TN.AdEx.C⁻ * TN.dt



        #===========================================
        						Learning Pfister
        		===========================================#
        if normalize
            if (tt % 200 == 0)
                norm = initial_synapses / sum(AMPAsynapses)
                AMPAsynapses .*= norm
            end
        end

        ## Duplet traces update before learning rule
        if spiked
            o1stdp += 1.0
        end
        for syn = 1:length(AMPAsynapses)
            r1stdp .+= exc_prespikes
            if exc_prespikes[syn]
                AMPAsynapses[syn] -= o1stdp * (triplet.A⁻₂ + triplet.A⁻₃ * r2stdp[syn])
            end
            if spiked
                AMPAsynapses[syn] += r1stdp[syn] * (triplet.A⁺₂ + triplet.A⁺₃ * o2stdp)
            end
            r2stdp .+= exc_prespikes
        end

        ## Triplet traces update after learning rule
        if spiked
            o2stdp += 1.0
        end

        r2stdp .*= exp(-TN.dt / triplet.τˣ)
        o2stdp *= exp(-TN.dt / triplet.τʸ)
        r1stdp .*= exp(-TN.dt / triplet.τ⁺)
        o1stdp *= exp(-TN.dt / triplet.τ⁻)
        AMPAsynapses[findall(AMPAsynapses .> max_efficacy)] .= max_efficacy
        AMPAsynapses[findall(AMPAsynapses .< 0.0)] .= 0

        #===========================================
        						Close learning
        		===========================================#
        ## close simulation step
        exc_prespikes .= false
        inh_prespikes .= false

        ## Recordings
        weights[:, tt] .= AMPAsynapses
        voltage[tt] = s.v
        θ[tt] = s.θ
        firing_rate[tt] = _firing_rate

        # ## Store results
        # rec_triplet[1,tt] = r1stdp[1]
        # rec_triplet[2,tt] = r2stdp[1]
        # rec_triplet[3,tt] = o1stdp
        # rec_triplet[4,tt] = o2stdp
        #############################
    end
    ## Current and Voltage plot
    return weights, voltage, θ, firing_rate
    # return voltage, rec_synapses, rec_triplet, weights
end

##
function simulate_vogels(
    s::TN.Soma,
    simTime::Int64,
    NMDAsynapses::Vector{Float32},
    AMPAsynapses::Vector{Float32},
    GABAsynapses::Vector{Float32},
    spiketimes::Vector{Int64},
    firing_neurons::Vector{Vector{Int}};
    max_efficacy = 5.0,
    normalize = false,
    kwargs...,
)

    ## Initialization and loop variables
    total_steps = round(Int, simTime / TN.dt)
    last_spike = -1
    spiked = (s.v > TN.AdEx.θ)

    ## Learning
    i_trace = zeros(Float32, size(GABAsynapses))
    rate = 0.0
    # e_trace = zeros(Float32, size(AMPAsynapses))
    # initial_synapses = sum(AMPAsynapses, dims=1)[1,:]
    # norm = initial_synapses

    ## Recording
    voltage = zeros(total_steps)
    firing_rate = zeros(total_steps)
    weights = zeros(size(GABAsynapses)..., total_steps)
    θ = zeros(total_steps)
    # rec_triplet =zeros(4,total_steps)

    ## Weights AMPA

    ## Input spikes

    reverse!(spiketimes)
    reverse!(firing_neurons)
    next_ft = pop!(spiketimes)
    next_neurons = pop!(firing_neurons)

    exc_prespikes = falses(size(AMPAsynapses))
    inh_prespikes = falses(size(GABAsynapses))
    inhibitory = Vector{Int}()
    excitatory = Vector{Int}()

    # for tt in ProgressBar(1:total_steps)
    @fastmath @inbounds for tt in ProgressBar(1:total_steps)

        if TN.update_AdEx!(s, spiked)
            spiked = true
            last_spike = tt
            s.θ = TN.postspike_vSTDP.A
        else
            # I_ps += - TN.dt * I_ps*TN.postspike_vSTDP.τz⁻
            s.θ += -TN.dt * (s.θ - TN.AdEx.θ) * TN.postspike_vSTDP.τA⁻
            # s.v  += TN.dt * TN.AdEx.C⁻ * I_ps
            spiked = false
        end


        ## Input
        while next_ft == tt
            excitatory = next_neurons[next_neurons.>0]
            inhibitory = abs.(next_neurons[next_neurons.<0])
            exc_prespikes[excitatory] .= true
            inh_prespikes[inhibitory] .= true

            for neuron in excitatory
                TN.exc_spike_plastic!(
                    s,
                    AMPA_eff = AMPAsynapses[neuron],
                    NMDA_eff = NMDAsynapses[neuron],
                )

            end
            for neurons in inhibitory
                TN.inh_spike!(s, eff = GABAsynapses[neuron])
            end
            ## Pull the queue
            if length(spiketimes) > 0
                next_ft = pop!(spiketimes)
                next_neurons = pop!(firing_neurons)
            else
                next_ft = -1
            end
        end

        #===========================================
        						Learning Vogels
        		===========================================#
        if normalize
            if (tt % 200 == 0)
                norm = initial_synapses / sum(AMPAsynapses)
                AMPAsynapses .*= norm
            end
        end

        ## Duplet traces update before learning rule
        rate -= dt * rate / iSTDP.τy
        i_trace .+= inh_prespikes .- dt * i_trace ./ iSTDP.τy

        if spiked
            rate += 1.0
        end
        for syn = 1:length(GABAsynapses)
            if inh_prespikes[syn]
                GABAsynapses[syn] += iSTDP.η * (rate - iSTDP.α)
            end
            if spiked
                GABAsynapses[syn] += iSTDP.η * i_trace[syn]
            end
        end

        AMPAsynapses[findall(AMPAsynapses .> max_efficacy)] .= max_efficacy
        AMPAsynapses[findall(AMPAsynapses .< 0.0)] .= 0

        GABAsynapses[findall(GABAsynapses .< iSTDP.j⁻)] .= iSTDP.j⁻
        GABAsynapses[findall(GABAsynapses .> iSTDP.j⁺)] .= iSTDP.j⁺

        #=====================
        		  	Close Learning
        		======================#


        ## close simulation step
        inh_prespikes .= false
        exc_prespikes .= false

        ## Recordings
        weights[:, tt] .= GABAsynapses
        voltage[tt] = s.v
        θ[tt] = s.θ
        firing_rate[tt] = rate

    end
    ## Current and Voltage plot
    return weights, voltage, θ, firing_rate
    # return voltage, rec_synapses, rec_triplet, weights
end
