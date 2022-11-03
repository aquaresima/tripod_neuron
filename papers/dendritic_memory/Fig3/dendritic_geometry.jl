"""
This file reproduces Fig. 1A.

In the file we compute the axial conductance of the dendritic compartments, for all the physiological condition (HUMAN, MOUSE) and all the dendritic lengths.
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"))


# Set the parameters for the grid search algorithm and compute βC
βC = (TN.AdExTripod.Er - TN.AdExTripod.θ)/TN.AdExTripod.θ
ds = 0.0:0.01:6
ls = 50:1.:550


# Obtain points for
function get_line(x)
    xx = reshape(x,ds.len, length(ls))
    z = map(x->findfirst(y-> y==1,x), eachcol(xx))
    return map(x->isnothing(x) ? ds[end] : ds[x],z)
end



Ri    = TN.HUMAN.Ri # same in MOUSE.Ri
βdistal_sup = Vector{Bool}()
βdistal_inf = Vector{Bool}()
for l_ in ls
    for d_ in ds
        ga = TN.G_axial(Ri=Ri,d=d_*TN.μm,l=l_*TN.μm)
        push!(βdistal_sup,ga.val>βC*TN.AdExTripod.gl)
        push!(βdistal_inf,ga.val>0.5βC*TN.AdExTripod.gl)
    end
end

βsup = get_line(βdistal_sup)
βinf = get_line(βdistal_inf)

θmouse = Vector{Bool}()
θhuman = Vector{Bool}()
Ri    = TN.MOUSE.Ri
Rd_h    = TN.HUMAN.Rd
Rd_m    = TN.MOUSE.Rd
βm = Vector{Bool}()
for l_ in ls
    for d_ in ds
        ga = TN.G_axial(Ri=Ri,d=d_*TN.μm,l=l_*TN.μm)
        gm = TN.G_mem(Rd=Rd_m,d=d_*TN.μm,l=l_*TN.μm)
        push!(θmouse,gm.val<ga.val)
        gm = TN.G_mem(Rd=Rd_h,d=d_*TN.μm,l=l_*TN.μm)
        push!(θhuman,gm.val<ga.val)
    end
end
θmouse= get_line(θmouse)
θhuman = get_line(θhuman)

##
p1 = plot(ylabel="Diameter (μm)", colorbar=false, xlabel="Length (μm)", ylims=(0,5.5), xlims=(50,550), frame=:box);
p1 = plot!(ls, βinf  , label="", lw=3 )
p1 = plot!(ls, βsup  , label="", lw=3 )
p1 = plot!(ls, θhuman, label="", lw=3 , ls=:dot,c=:grey)
p1 = plot!(ls, θmouse, label="", lw=3 , ls=:dot,c=:grey)
annotate!([(480,3.6,Plots.text(L"\frac{\beta}{2}<\frac{g_{ax}}{g^s_m}<\beta", :black,font))]);
annotate!([(480,5.,Plots.text(L"\frac{g_{ax}}{g^s_m}>\beta",:black))]);
annotate!([(480,1.4,Plots.text(L"\frac{g_{ax}}{g^s_m}<\frac{\beta}{2}",:black,font))]);
annotate!([(350, 2.,Plots.text("No spiking \nregime", :black,:top))]);
annotate!([(150, 4.2,Plots.text("Spiking\nregime", :black,:top))])
annotate!(250,0.95,Plots.text("Mouse \n"*L"g^d_{m}<g_{ax}",  :grey))
annotate!(250,4.9,Plots.text("Human \n"*L"g^d_{m}<g_{ax}", :grey))
savefig(p1,plotsdir("dendritic_memory","Fig3A.pdf"))