using DrWatson
@quickactivate "Tripod"

using TripodNeuron
using JLD2
using ProgressBars, Logging
using Plots

include(projectdir("scripts", "up_down", "updown_plots.jl"))
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
include(projectdir("scripts", "plots", "default_plots.jl"))
filename = "optimal_IE_μmem-55.0_soma_only.jld2"
file = datadir("up_down", "soma_only", filename)
data = load(file)

@unpack opt_kies, models, νs, min_AMPA, min_NMDA, min_AMPA_EQ, min_KUHN = data

#
plot()
for key in keys(opt_kies)
    kies = opt_kies[key]
    plot!(νs, kies[:, 1], clims = (0, 2), label = string(key))
end
p0 = plot!(
    xscale = :log,
    legend = :topleft,
    bg_color_legend = :transparent,
    ylabel = "I/E",
    xlabel = "Input rate (kHz)",
    title = "Balance E/I ratio",
)



plot()
for syn in syns
    filep = datadir("up_down", "soma_only", "full_activity_$syn.jld2")
    v = load(filep)["rate"][1, 1, 1, :]
    plot!(νs, v, label = string(syn))
end
p1 = plot!(
    xscale = :log10,
    legend = false,
    ylabel = "Firing rate",
    xlabel = "Input rate (kHz)",
)

plot()
for syn in syns
    filep = datadir("up_down", "soma_only", "full_activity_$syn.jld2")
    v = load(filep)["std_v"][1, 1, 1, :]
    plot!(νs, v, label = string(syn))
end
p2 = plot!(
    xscale = :log10,
    legend = false,
    ylabel = "STD membrane",
    xlabel = "Input rate (kHz)",
)


#
plotsNMDA = []
for syn in syns
    file = datadir("up_down", "soma_only", "critical_window_$syn.jld2")
    v = load(file)["data_model"]
    p = heatmap(
        TN.νs,
        βs,
        v.critical[1, :, :]',
        clims = (0, 50),
        c = :amp,
        colorbar = false,
    )
    push!(plotsNMDA, p)
    m = argmax(mean(v.critical[1, :, :], dims = 2))[1]
    plot!(title = string(syn))
    # vline!([TN.νs[m]])
    @show m
end
plot!(plotsNMDA[3], xlabel = "Input rate (kHz)", ylabel = "β")
plot!(plotsNMDA[4], xlabel = "Input rate (kHz)")
plot!(plotsNMDA[1], xticks = :none, ylabel = "β")
plot!(plotsNMDA[2], xticks = :none)
p3 = plot(
    plotsNMDA...,
    layout = (2, 2),
    size = (800, 600),
    legend = false,
    yticks = ([0, 250, 500]),
    bg_color_legend = :transparent,
    xscale = :log,
)
pp = colorbar(
    0,
    50,
    4,
    c = :amp,
    ylabel = "Critical window (mV)",
    size = (100, 400),
    orientation = :vertical,
    guidefontsize = 18,
)
layout = @layout [a{0.05w} _ b{0.9w}]
p3 = plot(pp, p3, layout = layout, size = (800, 600), margin = 2mm)
#




layout = @layout([
    a{0.3h}
    b{0.5w} c{0.5w}
    d{0.5h}
])
zz = plot(p0, p1, p2, p3, layout = layout, size = (800, 1000), lw = 4, margin = 5mm)

##
path = mkpath(plotsdir("up_down", "Figures"))
savefig(zz, joinpath(path, "FigApp2.pdf"))

plotsNMDA = []
for syn in syns
    file = datadir("up_down", "soma_only", "critical_window_$syn.jld2")
    v = load(file)["data_model"]
    p = heatmap(v.output_rates[:, :]', c = :amp, clims = (0, 30), colorbar = false)
    push!(plotsNMDA, p)
    m = argmax(mean(v.critical[1, :, :], dims = 2))[1]
    plot!(title = string(syn))
    # vline!([TN.νs[m]])
    @show m
end
plot(
    plotsNMDA...,
    layout = (2, 2),
    size = (800, 600),
    legend = false,
    yticks = ([0, 250, 500]),
)
