# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using DrWatson, Revise
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using Random


function input_protocol_tt(tt::Int64, comp::Int64 ;eff=1.f0, spikes::Matrix{Int})
    if spikes[comp,tt] >0
        return spikes[comp,tt]*eff
    else
        return 0
    end
end

##
Random.seed!(13)
simtime = 1000
model = TN.H_distal_proximal
label = TN.labels[3]
spikes = TN.make_balance(simtime,model=model, rate=5_000)
spikes[5:6,:] .*=0.70
AdExTripod = TN.AdExParams(Er=-70.6, u_r=-70.6)
voltage = TN.run_tripod(spikes, simtime; model..., do_spikes=true, AdEx=AdExTripod)
q= plot(voltage[1,:], c=:black, lw =3, label="Soma", title="Model: $label")
plot!(voltage[2,:], label="Dendrite: $(model.ds[1]) μm")
plot!(voltage[3,:], label="Dendrite: $(model.ds[2]) μm")
spike = TN.get_spike_times(voltage)[1]
plot!(xticks=([spike*10, 10*(spike+TN.AdEx.up)],[spike, spike+TN.AdEx.up]))
plot(q, xlims=((spike-60),(spike+100)))
annotate!(spike+60,0,Plots.text("BAP "*L"g_{ax}"*"= 1",  :black, :right))
p = plot!(size=(600,600), margin=3mm,)

annotate!(835, 45, Plots.text("A", 18, :black, :bold))
plot!(p, xlabel="", xticks=:none)

Random.seed!(13)
postspike=TN.PostSpike(BAP_gax=2,)
voltage = TN.run_tripod(spikes, simtime; model..., do_spikes=true, AdEx=AdExTripod,postspike=postspike)
q= plot(voltage[1,:], c=:black, lw =3, label="Soma",)
plot!(voltage[2,:], label="Dendrite: $(model.ds[1]) μm")
plot!(voltage[3,:], label="Dendrite: $(model.ds[2]) μm")
spike = TN.get_spike_times(voltage[1,:])[1]
plot!(xticks=([spike, (spike+TN.AdEx.up*10), spike+(TN.AdEx.up+TN.AdEx.idle)*10],floor.(Int,[spike, spike+TN.AdEx.up*10, spike+(TN.AdEx.up+TN.AdEx.idle)*10] ./10)), rotation=-45)
plot(q, xlims=((spike-60),(spike+100)))
# plot(q, xlims=(4200,4400))
q = plot!(size=(600,600), margin=3mm,)
plot!(xlabel="Time (ms)")
annotate!(spike+60,0,Plots.text("BAP "*L"g_{ax}"*"= 2",  :black, :right))

annotate!(835, 45, Plots.text("B", 18, :black, :bold))
p = plot(p,q, layout=(2,1), legend=:topleft, top_margin=10Plots.mm)

savefig(p,plotsdir("dendritic_memory","Revision_Backprop_comparison.pdf"))

##
Random.seed!(49)
simtime = 1000
model = TN.H_distal_proximal
label = TN.labels[3]
spikes = TN.make_balance(simtime,model=model, rate=5_000)
spikes[1,499] = 20.
AdExTripod = TN.AdExParams(Er=-70.6, u_r=-70.6)
voltage = TN.run_tripod(spikes, simtime; model..., do_spikes=true, AdEx=AdExTripod)
q= plot(voltage[1,:], c=:black, lw =3, label="Soma")#, title="Model: $label")
plot!(voltage[2,:], label="Dendrite: $(model.ds[1]) μm")
plot!(voltage[3,:], label="Dendrite: $(model.ds[2]) μm")
spike = TN.get_spike_times(voltage)[1]
plot!(xticks=([spike*10, 10*(spike+TN.AdEx.up)],[spike, spike+TN.AdEx.up]))
plot(q, xlims=((spike-60),(spike+100)))
p = plot!(size=(600,400), margin=3mm,)
plot!(xticks=([spike, (spike+TN.AdEx.up*10), spike+(TN.AdEx.up+TN.AdEx.idle)*10],floor.(Int,[spike, spike+TN.AdEx.up*10, spike+(TN.AdEx.up+TN.AdEx.idle)*10] ./10)), rotation=-45)

annotate!(835, 45, Plots.text("A", 18, :black, :bold))
ylims!((-100,40))
plot!(yticks=(-80:20:20, -80:20:20))
# plot!(p, xlabel="", xticks=:none)
plot!(xlabel="Time (ms)", ylabel="Membrane potential (mV)", legendfontsize=13)
savefig(p,plotsdir("dendritic_memory","Revision_Backprop.pdf"))