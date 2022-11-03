
"""
This file reproduces Fig. 1B.

In the file we compute the dendritic integration timescale of the dendritic compartments, for all the physiological condition (HUMAN, MOUSE) and all the dendritic lengths.
"""
## include Tripod
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"))

# Set the parameters for the grid search algorithm and compute βC
βC = (TN.AdExTripod.Er - TN.AdExTripod.θ)/TN.AdExTripod.θ
dr = 100:1:450
timescale = zeros(Float64,4,length(dr) )
colors    = zeros(RGBA,4,length(dr) )

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function get_timescale(pm)
    @unpack g_ax, Rm, C⁻ = pm
    C = 1/pm.τm⁻/pm.Rm
    return C/ (g_ax+1/Rm)
end

function get_color(params)
    @unpack g_ax = params
    G = βC*TN.AdExTripod.gl
    if 0.5G > g_ax
        return GREY
    elseif G > g_ax && 0.5G < g_ax
        return palette(:blues,2)[1]
    else
        return palette(:blues,2)[2]
    end
end

# get_color(params)
for n in eachindex(dr)
    for (a,S) in zip([1,2],["H","M"])
        d  = dr[n]
        pm   = TN.PassiveMembraneParameters("d2", S, 4, d)
        model ="$S. 1->s, $d,4"
        colors[a,n] = get_color(pm)
        timescale[a,n] = get_timescale(pm)
    end
    for (a,S) in zip([3,4],["H","M"])
        d  = dr[n]
        pm   = TN.PassiveMembraneParameters("d2", S, 2.5, d)
        colors[a,n] = get_color(pm)
        timescale[a,n] = get_timescale(pm)
    end
end
##
p = plot()
vline!([150], c=:black, ls=:solid, lw=5, alpha=0.3 ,label="");
vline!([400], c=:black, ls=:solid, lw=5, alpha=0.3 ,label="")
plot!(dr,timescale[1,:],c=colors[1,:],  lw=6, alpha=0.5,label="Human");
plot!(dr,timescale[2,:],c=colors[2,:],  lw=6, alpha=0.5,label="");
plot!(dr,timescale[3,:],c=colors[3,:],  lw=6, alpha=0.5,label="");
plot!(dr,timescale[4,:],c=colors[4,:],  lw=6, alpha=0.5,label="");
plot!(dr,timescale[2,:],color=:black, lw=3,ls=:dot,label="");
plot!(dr,timescale[4,:],color=:black, lw=3,ls=:dot,label="Mouse");
annotate!([(200,1.5, Plots.text("thin", 14, :black, :left))])
annotate!([(300,0.5, Plots.text("thick", 14, :black, :left))])
annotate!([(155,2.7, Plots.text("proximal", 14, :black, :left))])
annotate!([(335,2.7,Plots.text("distal", 14, :black, :left))])
plot!(xticks=([150,300,400],[150,300,400]), xlims=(100,500))
ylabel!("Membrane timescale  (ms)")
plot!(legendfontsize=12, bglegend=:transparent, 
        fglegend=:transparent, legend=:top)
xlabel!("dendritic length "*L"(\mu m)")
savefig(p,plotsdir("dendritic_memory","Fig3B.pdf"))