
##
using Distributed
println("Running Julia script with workers: $(workers())")
using DrWatson
@quickactivate "network_models"
using SpikingNeuralNetworks
SNN.@load_units;
using SNNUtils
using YAML
using Statistics
using Logging
using Plots
include("../models.jl")
##


## Set some test parameters. To test the stimulus is not very important
# begin 
#     @unpack models, kie, gsyn, min_AMPA, min_NMDA, νs = SNNUtils.tripod_balance.dend
#     syn_soma = TripodSomaSynapse
#     syn_dend = EyalEquivalentNAR(1.8)
# end

adex = AdExSoma( C=281pF, gl=40nS, Vr = -50.6, Er = -70.6, ΔT = 2, Vt = -50.4, a = 4, b = 80.5, τw = 144, up = 1ms, τabs= 1ms)
ps = SNN.PostSpike(A= 10.0,τA= 30.0)
ν = 3kHz
L = 300um
β = 200
nmda = true
Nd = 2
NAR = 0
Vs = -55mV


β =200
L = 250um
ν = 2.5kHz

E = Multipod(L, N=1, Nd=2,
    dend_syn = EyalEquivalentNAR(NAR),
    param= adex,
    postspike = ps
    )


param(β, ν, L) = SNN.BalancedStimulusParameter(
    kIE = SNNUtils.compute_kei(L, ν, NAR=NAR, Nd=Nd, Vs=Vs),
    wIE = 1.f0,
    r0 = ν,
    β =β,
)
stim = (
exc_stim_d1 = SNN.BalancedStimulus(E, :he, :hi, 1, cells=:ALL, μ=1.0, param = param(β, ν, L)),
exc_stim_d2 = SNN.BalancedStimulus(E, :he, :hi, 2, cells=:ALL, μ=1.0, param = param(β, ν, L)),
# exc_stim_d3 = SNN.BalancedStimulus(E, :he, :hi, 3, cells=:ALL, μ=1.0, param = param(β, n, l)),
)
model = SNN.merge_models(E=E, stim_d=stim )
SNN.sim!(model=model, duration = 2s,  pbar=true)
SNN.monitor(model.pop.E, [ :v_d, :v_s, :fire], sr=8000Hz)
SNN.sim!(model=model, duration = 20s,  pbar=true, dt=0.125)

#
st = spiketimes(model.pop.E)[1]
p =plot()
vecplot!(p,model.pop.E, :v_d, sym_id=1, r=0:2ms:20s, neurons=1)
vecplot!(p,model.pop.E, :v_d, sym_id=2, r=0:2ms:20s, neurons=1)
# vecplot!(p,model.pop.E, :v_d, sym_id=3, r=0:2ms:20s, neurons=1)
vecplot!(p,model.pop.E, :v_s, r=0:0.1ms:20s, neurons=1, c=:black)
scatter!(st/1000, 0.0*st .+ 5, c=:black, ms=4)
kIE = compute_kei(L, ν, NAR=NAR, Nd=Nd, Vs=Vs)
p1 = plot!(ylims=:auto, xlims=(14,16), title="Model: $Nd dendrites, $(round(Int, L*10000))um, NAR: $NAR. \nInput: ISICV $(β/1000+1), Exc rate: $(round(ν, digits=2))kHz, kIE: $(round(kIE, digits=2))")
# = istogram(getvariable(model.pop.E, :v_s)[1,:], bins=-70:1:-30, normed=true, c=:black, alpha=0.5)
p2= histogram(getvariable(model.pop.E, :v_s)[1,:], bins=-70:1:0, normed=true, c=:red, alpha=0.5, label="Soma")
p2= histogram!(getvariable(model.pop.E, :v_d)[1,2,:], bins=-70:1:0, normed=true, c=:blue, alpha=0.5, label="Dendrite")
plot(p1, p2, layout=(2,1), size=(800,800))
##

νs  = exp.(range(log(0.01), log(2), length=10))
plot()
for nmda in [true, false]
    for Nd in 1:3
        ν = 6kHz
        Vs = -55mV
        models = [150, 200, 250, 300, 350, 400] .*um
        vds = zeros(Float32, length(models))
        Threads.@threads for n in eachindex(models)
            NAR = nmda ? 1.8 : 0.0
            kie = compute_kei(models[n], ν, NAR=NAR, Nd=Nd, Vs=Vs)
            vs, vd = no_spike_model(l=models[n], ν=ν, Nd=Nd, nmda=nmda, kie = kie, β=0, dt=0.125ms/4)
            vds[n] = mean(vs)
        end
        name = nmda ? "NMDA_$Nd" : "AMPA_$Nd"
        plot!(models, vds, label=name, lw=3)
    end
end
plot!(ylims=(-80, -40), xlims=:auto)