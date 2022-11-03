using Plots, Revise, Random
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

function protocol(tt::Int64, comp::Int64;somaspike=false)
    if tt == 30
        comp == 1 && return 1.
        comp == 2 && return 1
        comp == 4 && return 1
        comp == 5 && return 1
    end
    return 0.
end

duration =150
default(palette=:tab10, lw=3)
synapses_human= TN.run_tripod((400,400), protocol, duration, rec_syn=true)
synapses_mouse= TN.run_tripod((100,100), protocol, duration, 
            syn_model=TN.mouse_synapses, rec_syn=true)

##
# synapses_dend_human
xx = TN.dt:TN.dt:duration
p=plot(xx, synapses_human[1,:,1], c=:black, label="AMPA", ls=:solid)
plot!(xx, synapses_mouse[2,:,2], c=RED, label="NMDA mouse Barrel Cortex", ls=:dash)
plot!(xx, synapses_human[2,:,2],c=BLU, label="NMDA Human L2/3", ls=:dot)
plot!(xlabel="Time (ms)", yaxis = false, lw=6)
q = plot(xx, synapses_human[2,:,4], c=:black, label="GABAb", ls=:solid)
plot!(xx, synapses_human[2,:,3], c=GREEN, label="GABAa dendrite", ls=:dash)
plot!(xx, synapses_human[1,:,3],c=BLU, label="GABAa soma", ls=:dot)
plot!(xlabel="Time (ms)", yaxis = false, lw=6)


plot(p,q, layout=(2,1), frame=:none)
savefig(plotsdir("dendritic_memory","Fig2A.pdf"))