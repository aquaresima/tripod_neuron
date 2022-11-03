"""
This code fits the parameters of Soma and Dendrite GABAa to the IPSP
reported in Miles 1996.

It recreates a stimuli consistent with the paper methods and use it to fit the IPSP.
"""
###
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"))
using HDF5, LsqFit

## Load data from figure
file = joinpath(@__DIR__,"fit_synapses","data_fit","miles_data_minimal.h5")
fid = h5open(file,"r")
miles_range = read(fid["x"])
soma_data = read(fid["soma"])
dendrite_data = read(fid["dend"])
close(fid)

##
function protocol(tt::Int64, comp::Int64; target_comp)
    if tt ==SPIKETIME && comp == target_comp
        return SYN_CONTACT/2
    end
    return 0
end

SYN_CONTACT = 10 # Number of synapses: 5 per soma, and 5 per dendrites 
Er = TN.AdEx.Er
xrange= 0.1:0.1:127
SPIKETIME = 1270 # 27ms

MilesTestDend = TN.GABAergic(
 TN.Receptor(E_rev = Er-30,τr= 8, τd= 29., g0= 0.27),
 TN.Receptor(E_rev = -120,τr= 30, τd= 400., g0=0.006)
)
MilesTestSoma = TN.GABAergic(
 TN.Receptor(E_rev = Er-30,τr= 0.1, τd= 15., g0= 0.50),
 TN.Receptor()
)
fit_synapses = TN.TripodSynapses(
	Esyn_soma = TN.Synapse(TN.DuarteGluSoma, MilesTestSoma),
	Esyn_dend = TN.Synapse(TN.EyalGluDend, MilesTestDend),
)
model = (ds = (200,200), syn_model = fit_synapses , species="H")
voltage = TN.run_tripod(model.ds, protocol,227; model..., target_comp=5)
s = plot(xrange,.+ voltage[1,1001:end] .- voltage[1,1], color=RED, label="",lw=3)
plot!(miles_range, dendrite_data, color=RED, linestyle=:dash, label="",lw=3)

model = (ds = (300,300), syn_model = fit_synapses , species="H")
voltage = TN.run_tripod(model.ds, protocol,227; model..., target_comp=4)
s = plot!(xrange,.+ voltage[1,1001:end] .- voltage[1,1], color=BLU, label="",lw=3)
plot!(miles_range, soma_data, color=BLU,linestyle=:dash, label="",lw=3)

plot!(s, [],[],color=:black, linestyle=:dash, label = "Miles et al. 1996");
plot!(s, [],[],color=RED, linestyle=:solid, label = "Dendritic spike")
plot!(s, [],[],color=BLU, label = "Perisomatic spike");
s = plot!(s,[50,70], [.25, .25], color=:black, label=false);
s = plot!(s,[130,130], [0.25, -0.75], color=:black, label=false)
annotate!([(133,-0.25, Plots.text("1 mV",:left, 18)),(60,.40,Plots.text("20 ms",19))])
plot!(xlims=(0,155), ylims=(-1.5, 0.9))
# annotate!([(1700,-40, Plots.text("Dendrite")),(1700,-50,Plots.text("Soma"))]);
# plot!(s,xlims=(1500,2600), ylims=(-55,-48))
plot!(s,xaxis=false, yaxis=false, legend=:topright, legendfontsize=13, ticks=:none);
plot(s, linewidth=3.)

savefig(s, plotsdir("dendritic_memory","Fig2B.pdf"))

s
## Fit parameters
