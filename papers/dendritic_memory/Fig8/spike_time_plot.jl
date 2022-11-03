using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using StatsBase, Statistics, JLD 
using MLDataUtils, MLJLinearModels, RollingFunctions
using NamedTupleTools
file =datadir("dendritic_memory","spike_time_data.jld")

##
data = JLD.load(file) |> adict ->namedtuple(Symbol.(keys(adict)), values(adict))
@unpack dend, soma, no_inh, base_values = data
reverse_amp=[]
for a in 1:256
    x = reverse(ColorSchemes.amp)[a]
    push!(reverse_amp,RGBA(x,1))
end
reverse_amp = cgrad(reverse_amp);

heatmap(mean(no_inh, dims=1)[1,:,1,:], clims=(5,40), c=reverse_amp, xlabel="Retrieval delay (s)", ylabel="Excitation (kHz)", title="Memory access", size=(400,400),
    yticks=(range(1,26,3), range(0,5,3)),
    xticks=(range(1,26,3), range(0,5,3)),
)


p = plot(
    heatmap(mean(dend,dims=1)[1,:,:,1],c=reverse_amp, clims=(0,30)),
    heatmap(mean(dend,dims=1)[1,:,:,2],c=reverse_amp, clims=(0,30)),
    heatmap(mean(dend,dims=1)[1,:,:,3],c=reverse_amp, clims=(0,30)),
    heatmap(mean(soma,dims=1)[1,:,:,1],c=reverse_amp, clims=(0,30),
    xlabel="Inhibition (kHz)", ylabel="Excitation (kHz)"),
    heatmap(mean(soma,dims=1)[1,:,:,2],c=reverse_amp, clims=(0,30)),
    heatmap(mean(soma,dims=1)[1,:,:,3],c=reverse_amp, clims=(0,30)
    ),
    colorbar=false,
    yticks=(range(1,26,3), range(0,5,3)),
    xticks=(range(1,26,3), range(0,5,3)),
    size=(800,600),

    margin=4mm
)
savefig(p,plotsdir("dendritic_memory","Fig8C.pdf"))


mean(dend,dims=1)
p = plot(
    heatmap(mean(dend,dims=1)[1,:,:,1]-mean(soma,dims=1)[1,:,:,1],c=:redblue, clims=(-20,20)),
    heatmap(mean(dend,dims=1)[1,:,:,2]-mean(soma,dims=1)[1,:,:,2],c=:redblue, clims=(-20,20)),
    heatmap(mean(dend,dims=1)[1,:,:,3]-mean(soma,dims=1)[1,:,:,3],c=:redblue, clims=(-20,20), xlabel="Inhibition (kHz)", ylabel="Excitation(kHz)"),
    layout=(3,1), 
    yticks=(range(1,26,3), range(0,5,3)),
    xticks=(range(1,26,3), range(0,5,3)),
    size=(600,900),
    margin=4mm
)
savefig(p,plotsdir("dendritic_memory","Fig8D.pdf"))