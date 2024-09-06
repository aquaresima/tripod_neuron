"""
This file reproduces Fig 2A and 2B of the paper
It is necessary to run compute_dendritic_length.jl before it
"""

using DrWatson
@quickactivate "Tripod"
using Revise
using Plots, JLD, HDF5, Statistics, EllipsisNotation, RollingFunctions
using TripodNeuron
include(projectdir("scripts", "plots", "default_plots.jl"));

## Read data
data = load(datadir("up_down", "rates", "dendritic_length_ds.jld2"))

rs = data["data"][1, ..]
cv = data["data"][2, ..]
membrane = data["data"][3:end, ..]
νs = data["νs"]
ds = data["ds"]
data
# @info "Loaded data inh_ratio: $(data["inh_ratio"]), Ereset: $(data["reset"])"

plots_rs = []
plots_cv = []
for x in [15, 25, 30, 35]
    @show x

    p = heatmap(
        ds[2:end],
        ds[2:end],
        rs[2:end, 2:end, x],
        yticks = :none,
        xticks = :none,
        frame = :box,
        colorbar = false,
        title = L"ν_{exc}" * " = $(round(νs[x], digits=1)) kHz",
        titlefontsize = 18,
        clims = (0, 70),
        c = :roma,
    )
    push!(plots_rs, p)
    q = heatmap(
        ds[2:end],
        ds[2:end],
        (cv[2:end, 2:end, x]),
        yticks = :none,
        xticks = :none,
        frame = :box,
        colorbar = false,
        title = L"ν_{exc}" * " = $(round(νs[x], digits=1)) kHz",
        titlefontsize = 15,
        clims = (0, 2.5),
        c = :lajolla,
    )
    push!(plots_cv, q)
    # plot!(p,xticks=((100,500),(100,500)))
    # plot!(q,xticks=((100,500),(100,500)))
end

_xticks = ([100, 300, 500], [100, 300, 500])
plot!(
    plots_cv[1],
    xlabel = "Dendritic length " * L"(\mu" * "m)",
    ylabel = "Dendritic length " * L"(\mu" * "m)",
)
plot!(plots_cv[1], xticks = _xticks, yticks = _xticks)
plot!(plots_rs[1], ylabel = " ", xticks = _xticks, yticks = _xticks)
# layout = @layout [ a{0.05w} grid(1,4)]
layout = @layout [grid(1, 4) a{0.05w}]
hm = colorbar(
    0,
    70;
    c = :roma,
    title = "Firing rate (Hz)",
    topmargin = 20Plots.mm,
    bottommargin = 20Plots.mm,
    leftmargin = 25Plots.mm,
)
annotate!(hm, -2.5, 120, text("A", :center, 20))
p = plot(plots_rs..., hm, layout = layout, titlefontsize = 15)
hm = colorbar(
    0,
    2.5,
    3;
    c = :lajolla,
    title = "ISI CV (Hz)",
    topmargin = 20Plots.mm,
    bottommargin = 20Plots.mm,
    leftmargin = 20Plots.mm,
    digits = 2,
)
layout = @layout [grid(1, 4) a{0.05w}]
annotate!(hm, -2.8, 150, text("B", :center, 20))
q = plot(plots_cv..., hm, layout = layout)
layout = @layout [
    _
    a{0.45h}
    b{0.45h}
    _
    _
]
zz1 = plot(p, q, layout = layout, size = two_columns_size, margin = 5Plots.mm)
zz1
