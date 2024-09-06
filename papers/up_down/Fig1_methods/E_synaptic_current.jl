
using DrWatson
@quickactivate "Tripod"

using TripodNeuron
using Statistics, StatsPlots, JLD2, Plots, RollingFunctions
include(projectdir("scripts", "plots", "default_plots.jl"))


## This plot shows that the model is not bistable (Fig 1B)
##
syn = TN.human_synapses.Esyn_dend

NARSyn(NAR) = Synapse(NARGluDend(NAR), MilesGabaDend)

function current(v; NAR = 1.8, KIR = false)
    if KIR
        return TN.syn_current_GABABKIR(v, ones(Float32, 4), NARSyn(NAR))
    else
        return TN.syn_current(v, ones(Float32, 4), NARSyn(NAR))
    end
end

default(linewidth = 3)
plot()
hline!([0], c = :gray, label = "")
plot!(
    -90:0,
    current.(-90.0f0:1.0f0:0.0f0),
    c = :black,
    label = "Tripod model (Quaresima 2022)",
)
plot!(
    -90:0,
    current.(-90.0f0:1.0f0:0.0f0, KIR = true),
    c = :black,
    label = "GABAb/KIR model (Sanders 2013) ",
    ls = :dash,
)
nar_values = 1.8:2:20
colors = range(mpi_palette[1], mpi_palette[2], length = length(nar_values))
for (nar, c) in zip(nar_values, colors)
    plot!(-90:0, current.(-90.0f0:1.0f0:0.0f0, NAR = nar), c = c, ls = :dash, label = "")
end
p = plot!(-90:0, current.(-90.0f0:1.0f0:0.0f0), c = :black, label = "")

# heatmap(p,zeros(10,10), subplot=2, inset=inset)

plot!(
    p,
    legend = :topright,
    ylabel = "Synaptic current (nA)",
    xlabel = "Membrane potential (mV)",
    frame = :box,
    ylims = (-300, 100),
    bglegend = :transparent,
)

inset = bbox(0.22, 0.65, 0.2, 0.05)
p = colorbar!(
    p,
    inset,
    1.8,
    20,
    c = range(mpi_palette[1], mpi_palette[2], length = length(nar_values)),
    digits = 1,
    horizontal = true,
    topmargin = 5Plots.mm,
    yrotation = 45,
    xlabel = "NAR (nA)",
    guidefontsize = 12,
    tickfontsize = 15,
    labelfontsize = 18,
)

##

savefig(p, joinpath(updown_plot, "E_synaptic_current.pdf"))
