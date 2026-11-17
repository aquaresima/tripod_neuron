# Clopath 2010
struct STDP
    #voltage based stdp
    a⁻::Float32    #ltd strength (pF/mV) # a*(V-θ) = weight
    a⁺::Float32    #ltp strength (pF/mV)
    θ⁻::Float32 #ltd voltage threshold (mV)
    θ⁺::Float32 #ltp voltage threshold (mV)
    τu⁻::Float32  #timescale for u variable   (1/ms)
    τv⁻::Float32  #timescale for v variable   (1/ms)
    τx⁻::Float32  #timescale for x variable   (1/ms)
    ϵ::Float32  # filter for delayed membrane potential.
    j⁻::Float32  # minimum weight
    j⁺::Float32  # maximum weight

end

#Vogel 2011
#inhibitory stdp
struct ISTDP
    τy::Float32
    η::Float32
    r0::Float32
    α::Float32
    j⁻::Float32  # minimum weight
    j⁺::Float32  # maximum weight
end

struct TripletRule
    A⁺₂::Float32
    A⁺₃::Float32
    A⁻₂::Float32
    A⁻₃::Float32
    τˣ::Float32
    τʸ::Float32
    τ⁺::Float32
    τ⁻::Float32
end

# Gutig 2003
struct NLTAH
    τ::Float32
    λ::Float32
    μ::Float32
end


# const postspike_10 = PostSpike(10, round(Int,1/dt), 1/30)

struct LearnParams
    inh::Bool
    exc::Bool
    function LearnParams(; inh::Bool, exc::Bool)
        new(inh, exc)
    end
end


function lkd_stdp()
    altd = 0.0008 #ltd strength (pF/mV) # a*(V-θ) = weight
    altp = 0.0014 #ltp strength (pF/mV)
    thetaltd = -70.0 #ltd voltage threshold (mV)
    thetaltp = -49.0 #ltp voltage threshold (mV)

    tauu = 10.0 #timescale for u variable   (ms)
    tauv = 7.0  #timescale for v variable   (ms)
    taux = 15.0 #timescale for x variable   (ms)
    ϵ = 0.0   # filter for delayed voltage

    jeemin = 0.78 #minimum ee strength (pF)  pF/ms = nS #xerise/difEtau = ge
    jeemax = 11.4 #maximum ee strength (pF)
    return STDP(
        altd,
        altp,
        thetaltd,
        thetaltp,
        1 / tauu,
        1 / tauv,
        1 / taux,
        ϵ,
        jeemin,
        jeemax,
    )
end

function clopath_vstdp_visualcortex()
    altd = 14e-5 / 2.5 #ltd strength (pF/mV) # a*(V-θ) = weight
    altp = 8e-5 / 2.5 #ltp strength (pF/mV)
    thetaltd = -70.6 #ltd voltage threshold (mV)
    thetaltp = 25.3 #ltp voltage threshold (mV)

    tauu = 10.0 #timescale for u variable   (ms) τ⁻
    tauv = 7.0  #timescale for v variable   (ms) τ⁺
    taux = 15.0 #timescale for x variable   (ms)
    ϵ = 1.0   # filter for delayed voltage

    jeemin = 0.5 #minimum ee strength (pF)  pF/ms = nS #xerise/difEtau = ge
    jeemax = 1.5 #maximum ee strength (pF)

    return STDP(
        altd,
        altp,
        thetaltd,
        thetaltp,
        1 / tauu,
        1 / tauv,
        1 / taux,
        ϵ,
        jeemin,
        jeemax,
    )

end

function bono_vstdp()
    altd = 5e-4 / 2.5 #ltd strength (pF/mV) # a*(V-θ) = weight
    altp = 0.0015 / 2.5 #ltp strength (pF/mV)
    thetaltd = -69.0 #ltd voltage threshold (mV)
    thetaltp = -15.0 #ltp voltage threshold (mV)

    tauu = 35.0 #timescale for u variable   (ms) τ⁻
    tauv = 35.0  #timescale for v variable   (ms) τ⁺
    taux = 15.0 #timescale for x variable   (ms)
    ϵ = 1.0   # filter for delayed voltage

    jeemin = 0.5 #minimum ee strength (pF)  pF/ms = nS #xerise/difEtau = ge
    jeemax = 1.5 #maximum ee strength (pF)

    return STDP(
        altd,
        altp,
        thetaltd,
        thetaltp,
        1 / tauu,
        1 / tauv,
        1 / taux,
        ϵ,
        jeemin,
        jeemax,
    )

end

function lkd_istdp()
    ## Inhibition
    tauy = 20.0 #decay of inhibitory rate trace (ms)
    eta = 1.0   #istdp learning rate    (pF⋅ms) eta*rate = weights
    r0 = 0.003   #target rate (khz)
    alpha = 2 * r0 * tauy #rate trace threshold for istdp sign (kHz) (so the 2 has a unit)
    jeimin = 48.7 #minimum ei strength (pF)
    jeimax = 243 #maximum ei strength   (pF)

    return ISTDP(tauy, eta, r0, alpha, jeimin, jeimax)
end

function vogels_istdp()
    ## Inhibition
    tauy = 20.0 #decay of inhibitory rate trace (ms)
    eta = 1e-1   #istdp learning rate    (pF⋅ms) eta*rate = weights
    r0 = 0.005   #target rate (khz)
    alpha = 2 * r0 * tauy #rate trace threshold for istdp sign (kHz) (so the 2 has a unit)
    jeimin = 0.1 #48.7 #minimum ei strength (pF)
    jeimax = 243 #maximum ei strength   (pF)

    return ISTDP(tauy, eta, r0, alpha, jeimin, jeimax)
end

function pfister_visualcortex(alltoall::Bool = true, full::Bool = true)
    if alltoall
        if full
            return TripletRule(5e-10, 6.2e-3, 7e-3, 2.3e-4, 101.0, 125.0, 16.8, 33.7)
        else
            return TripletRule(0.0, 6.5e-3, 7.1e-3, 0.0, -1.0, 125.0, 16.8, 33.7)
        end
    else
        if full
            return TripletRule(8.8e-11, 5.3e-2, 6.6e-3, 3.1e-3, 714.0, 40.0, 16.8, 33.7)
        else
            return TripletRule(0.0, 5.2e-2, 8.e-3, 0.0, -1.0, 40.0, 16.8, 33.7)
        end
    end
end

function clopath_vstdp_neuron()
    ## membrane params
    C = 281        #(pF)
    gL = 30         #(nS) leak conductance #BretteGerstner2005 says 30 nS
    τm = C / gL       #(ms) membrane RC time

    # thresholds
    u_reset = -60.6      #(mV)
    Erest = -70.4
    θ = -50.4

    ## adaptation
    τw = 144        #(ms) adaptation time constant (~Ca-activated K current inactivation)
    a = 4.0        #(nS) 'sub-threshold' adaptation conductance
    b = 0.805       #(pA) 'sra' current increment

    up = 1.0
    idle = 2.0

    return AdExParams(τm, gL, Erest, u_reset, θ, τw, a, b, up, idle)
end

function bonoclopath_vstdp_neuron()
    ## membrane params
    C = 281        #(pF)
    gL = 40         #(nS) leak conductance #BretteGerstner2005 says 30 nS
    τm = C / gL       #(ms) membrane RC time

    # thresholds
    u_reset = -55.0      #(mV)
    Erest = -69.0
    θ = -50.4

    ## adaptation
    τw = 144        #(ms) adaptation time constant (~Ca-activated K current inactivation)
    a = 0.0        #(nS) 'sub-threshold' adaptation conductance
    b = 0.0       #(pA) 'sra' current increment

    return AdExParams(τm, gL, Erest, u_reset, θ, τw, a, b, up, idle)
end
