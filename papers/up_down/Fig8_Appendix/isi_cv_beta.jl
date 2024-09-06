##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
using EllipsisNotation

include(projectdir("scripts", "plots", "default_plots.jl"))
pyplot()
gr()

##
get_data(title, type = "data") =
    load(datadir("up_down", "bistability", "critical_beta_tau_$(title)_rate=12.jld2"))
NMDA = get_data("NMDA")["data"]
τs = get_data("NMDA")["τs"]
βs = get_data("NMDA")["βs"]
ν = get_data("NMDA")["ν"]
@info "Loaded data, $(keys(collect(values(NMDA))[1]))"
titles = ["distal\ndistal", "distal\nproximal", "proximal\nproximal", "soma\nonly"]
##

@unpack νs, ds = TN
ampa_color = :darkred
nmda_color = :darkblue
##


data = collect(values(load(datadir("up_down", "inputs", "cv_inputs.bson"))))[1]
@unpack cvs, τs, βs = data

q2 = contour(
    βs[2:end],
    τs,
    cvs[:, 2:end, 1],
    title = "CV of input spike trains",
    xlabel = "Fluctuations (β)",
    ylabel = "τ (ms)",
    legend = :topleft,
    yscale = :log10,
    xscale = :log10,
    fill = true,
    c = :amp,
    rightmargin = 10mm,
    lw = 0.5,
)
hline!([50], c = :black, lw = 2, ls = :dash, label = "τ=50ms")
annotate!([(50, 50, text("τ=50ms", 15, :bottom, :black))])
# hline!([1], c=:black, lw=2, ls=:dash,label="τ=50ms")
plotsNMDA = []
for (n, (k, v)) in enumerate(NMDA)
    p = contour(
        βs[2:end],
        τs,
        v.critical[1, 1:end, 2:end],
        clims = (0, 50),
        c = :amp,
        fill = true,
        lw = 0.5,
    )
    push!(plotsNMDA, p)
    m = argmax(mean(v.critical[1, :, :], dims = 2))[1]
    plot!(title = "Bimodality in distal-proximal")
    # vline!([TN.νs[m]])
    hline!([50], c = :black, lw = 2, ls = :dash, label = "τ=50ms")
    # hline!([1], c=:black, lw=2, ls=:dash,label="τ=50ms")
    annotate!([(50, 50, text("τ=50ms", 15, :bottom, :black))])
    @show m
end
p1 = plot(
    q2,
    plotsNMDA[2],
    yscale = :log10,
    xscale = :log10,
    layout = (1, 2),
    size = (900, 400),
    xlabel = "Fluctuations (β)",
    ylabel = "τ (ms)",
    legend = false,
    bottommargin = 10Plots.mm,
)

# plot(mean(NMDA["proximal-proximal"].critical[1,:,1:20], dims=2))
# τs

##

using Statistics
stt = 2
c = collect(palette(:viridis, 51 - stt))
# plotly()
# q3 = scatter(repeat(βs[2:end],inner=49), cvs[stt:end,2:end,1][:],NMDA["proximal-proximal"].critical[1,stt:end,2:end][:], c=c, label="", legend=false, xlabel="β ", zlabel="Bimodality index", ylabel="ISI CV", title="Bimodality vs CV", size=(400,400))
# plot!(guidefontsize=11, tickfontsize=9)
q3 = scatter(
    cvs[stt:end, 2:end, 1][:],
    NMDA["proximal-proximal"].critical[1, stt:end, 2:end][:],
    c = c,
    msc = c,
    label = "",
    legend = false,
    ylabel = "Bimodality index",
    xlabel = "ISI CV",
    title = "Bimodality vs CV",
    size = (400, 400),
)
cc = round(
    corspearman(
        NMDA["distal-distal"].critical[1, stt:end, 2:end][:],
        cvs[stt:end, 2:end, 1][:],
    ),
    digits = 2,
)
annotate!(2.5, 10, text("ρ=$cc", 15, :bottom, :black))
xs = round.(collect(exp.(range(log(τs[stt]), log(1000), length = 5))), digits = 2)
q4 = colorbar(
    values = xs,
    c = :viridis,
    ylabel = "τs (ms)",
    legend = false,
    size = (400, 400),
    guidefontsize = 16,
)
layout = @layout [a{0.7w} _ b{0.2w}]
p2 = plot(q3, q4, layout = layout, size = (600, 400), guidefontsize = 16)

zz = plot(
    p1,
    p2,
    layout = (2, 1),
    size = (900, 800),
    legend = false,
    bottommargin = 10Plots.mm,
)

savefig(zz, plotsdir("up_down", "Figures", "Fig9_Appendix.pdf"))
