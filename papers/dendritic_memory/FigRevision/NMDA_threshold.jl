using DrWatson, Revise
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

## values for paper table:
nmda(x) = TN.NMDA_nonlinear(TN.human_synapses.Esyn_dend.NMDA,x)*TN.human_synapses.Esyn_dend.NMDA.g0
function get_timescale(pm)
    @unpack g_ax, Rm, C⁻ = pm
    C = 1/pm.τm⁻/pm.Rm
    return C/ (g_ax+1/Rm)
end


for (a,S) in zip([1,2],["H","M"])
    for d in [400,150]
    pm   = TN.PassiveMembraneParameters("d2", S, 4, d)
    t = round(get_timescale(pm), digits=2)
    print(t, S, " & ")
    end
end
default(lw=3)
function my_min_nmda(target)
    μ = 1
    z = 1
    for x in -90f0:0.1f0:0f0
        if abs(nmda(x) -target) < μ
            μ =  abs(nmda(x) -target )
            z = x
        end
        # @show μ, x
    end
    return z
end
ampa_g = TN.human_synapses.Esyn_dend.AMPA.g0
nmda_g = TN.human_synapses.Esyn_dend.NMDA.g0
x0 = TN.AdExTripod.θ
x1 = my_min_nmda(ampa_g/3.3)
x2 = my_min_nmda(ampa_g)
plot()
plot!(-95f0:0.1f0:40f0,nmda.(-95f0:0.1f0:40f0), c=:black,label="NMDAR conductance")
hline!([nmda_g], c=:black, ls = :dot,lw=3, label="NMDAR gsyn")
plot!([-90,x1],[ampa_g/3.3, ampa_g/3.3], ls=:dot,lw=3, label="30% AMPAR gsyn", c=GREEN)
plot!([-90,x2],[ampa_g, ampa_g], c=BLU, ls=:dot, lw=3, label="AMPAR gsyn")
plot!([x1,x1],[0,nmda(x1)], c=GREEN, ls=:dot, label="", lw=3)
plot!([x2,x2],[0,nmda(x2)], c=BLU, ls=:dot, label="", lw=3)
annotate!(-55,1.4,Plots.text("NMDAR "*L"\bar gsyn",  :black, :right))
annotate!(-55,0.83,Plots.text("AMPAR "*L"\bar gsyn",  :black, :right,-90))
annotate!(-55,0.32,Plots.text(L"\frac{1}{3}"*"AMPAR "*L"\bar gsyn",  :black, :right))
vline!([TN.AdExTripod.θ], c=RED, ls=:dash)
annotate!(-50,1.6,Plots.text("AdEx spike\n threshold",  RED, :bottom))
vline!([0.f0], c=GREY, ls=:dash)
annotate!(-0,1.6,Plots.text("GluRs reverse \npotential",  GREY, :bottom))
plot!(-90f0:0.1f0:0f0,nmda.(-90f0:0.1f0:0f0), c=:black,label="")
plot!(legend=:topleft, lw=4, xlabel="Membrane potential (mV)", ylabel="NMDAR conductance (nS)")
annotate!(-19,1.,Plots.text("NMDAR",  :black, rotation = 0))
plot!(xticks=([-80, x0,x1,x2,0,20], [-80, x0,x1,x2,0,20]), xlims=(-90,25), rotation=45, legend=false, top_margin=15Plots.mm)

annotate!(-110, 1.7, Plots.text("A", 18, :black, :bold))
nmda(x) = round(TN.NMDA_nonlinear(TN.human_synapses.Esyn_dend.NMDA,x)*TN.human_synapses.Esyn_dend.NMDA.g0,digits=2)
p = plot!(yticks=([ nmda(x1),nmda(x2),nmda_g], [nmda(x1),nmda(x2),nmda_g]), ylims=(0,1.54))
plot!(size=(600,400))
savefig(p,plotsdir("dendritic_memory","Revision_NMDA_threshold.pdf"))
p