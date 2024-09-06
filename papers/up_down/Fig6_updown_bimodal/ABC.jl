
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Plots
using JLD2
using ProgressBars, Logging, EllipsisNotation
include(projectdir("scripts", "up_down", "updown_plots.jl"))
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
include(projectdir("scripts", "plots", "default_plots.jl"))
##

default(titlefontsize = 13)
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
simtime = 25_000
rate = 13;
β = 200;
models_hist = []
nmda = true
ratio = 1
for (model, title) in zip(TN.models, TN.labels)
    cond = TN.get_balance_conditions(model.ds..., rate, nmda = nmda)

    inputs = TN.make_spikes(simtime, β; cond...)
    voltage =
        TN.run_tripod(inputs, simtime; model..., AdEx = AdEx, syn_model = cond.syn_model)

    kde = globalKDE(2, voltage[1, :], v_range = -70:40)
    m = maximum(kde)
    c = critical_window(voltage[1, :])
    aph = APH(voltage, "", 2, yticks = :none)
    plot!(aph, ylims = (0, m + m / 4))
    push!(models_hist, aph)
    annotate!(-70, m + m / 3, text(title, :left, 11))
end
for p in models_hist[1:end-1]
    plot!(p, xticks = false)
end
plot!(models_hist[1], top_margin = 5Plots.mm)
plot!(
    models_hist[end],
    xlabel = "                                    Membrane potential (mV)",
    xrotation = 45,
)
q = plot(
    models_hist...,
    layout = (4, 1),
    yaxis = false,
    tickfontsize = 13,
    guidefontsize = 18,
    frame = :axes,
    yticks = :none,
    size = (400, 600),
    margin = 2Plots.mm,
)


model = TN.models[2]
cond = TN.get_balance_conditions(model.ds..., rate, nmda = nmda)
inputs = TN.make_spikes(simtime, β; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., AdEx = AdEx, syn_model = cond.syn_model)

p = Plots.histogram(
    voltage[1, :],
    normalize = true,
    bins = range(-80, stop = -40, length = 100),
    color = :black,
    label = "",
    title = "Kernel Density Estimate",
    titlefontsize = 13,
    ylims = (0, 0.5),
)
# aph = APH(voltage,"",3, yticks=:none)

# w_plots = []
c = palette(:reds, 6)
for (n, w) in enumerate([2, 4, 8, 16, 32, 50])
    kde = globalKDE(w, voltage[1, :], v_range = -80:-40)
    plot!(
        -80:-40,
        kde,
        lw = 2;
        c = c[n],
        label = "$w   mV",
        legendtitle = "Window",
        legend_fontsize = 19,
        legend = :topleft,
    )
end

m = maximum(globalKDE(1, voltage[1, :], v_range = -80:-40))
plot!(
    yticks = :none,
    frame = :axes,
    yaxis = false,
    margin = 2Plots.mm,
    ylims = (0, m * 1.5),
    xrotation = 45,
)

hm = zeros(50, 1)
hm[:, 1] .= collect(1:length(hm))
inset = (1, bbox(-0.6, 0.60, 0.1, 0.3, :bottom, :right))
xs = round.(Int, [1, 25, 50])
_xticks = (xs, xs)
q1 = heatmap!(
    hm,
    c = :amp,
    frame = :box,
    xticks = :none,
    cbar = false,
    inset = inset,
    subplots = 2,
    yticks = _xticks,
    ylabel = "Window",
    guidefontsize = 11,
    tickfontsize = 11,
)
z1 = plot(q, p, size = (600, 600))

##

get_data(syn) = load(datadir("up_down", "bistability", "critical_window_$syn.jld2"))

AMPA = get_data("AMPA")["data"]
NMDA = get_data("NMDA")["data"]
βs = get_data("NMDA")["βs"]

titles = ["distal\ndistal", "distal\nproximal", "proximal\nproximal", "soma\nonly"]

# keys(collect(values(NMDA))[1])



plotsNMDA = []
for (n, (k, v)) in enumerate(NMDA)
    p = heatmap(TN.νs, βs, v.critical[2, :, :]', clims = (0, 50), c = :amp)
    push!(plotsNMDA, p)
    m = argmax(mean(v.critical[2, :, :], dims = 2))[1]
    plot!(title = titles[n])
    # vline!([TN.νs[m]])
    @show m
end

plotsAMPA = []
for (n, (k, v)) in enumerate(AMPA)
    p = heatmap(TN.νs, βs, v.critical[2, :, :]', clims = (0, 50), c = :amp)
    push!(plotsAMPA, p)
    m = argmax(mean(v.critical[2, :, :], dims = 2))[1]
    plot!(title = titles[n])
    # vline!([TN.νs[m]])
    @show m
end
[plot!(plotsNMDA[x], xticks = :none, yticks = :none) for x in [1, 2, 4]]
[plot!(plotsAMPA[x], xticks = :none, yticks = :none) for x in [1, 2, 4]]
plot!(
    plotsNMDA[3],
    ylabel = "\nFlutctuation size (β)",
    xlabel = "                  Input rates (kHz)",
)
plot!(plotsAMPA[3], ylabel = " ", xlabel = " ", yticks = :none)
p = plot(plotsNMDA..., xscale = :log, clims = (0, 50))
q = plot(plotsAMPA..., xscale = :log, clims = (0, 50))
z2 = plot(p, q, size = (800, 600), colorbar = false, margin = 2Plots.mm, legend = false)
layout = @layout [a{0.95h,0.35w} _ a{0.6w} _]
zz1 = plot(
    z1,
    z2,
    size = two_columns_size .* (1.2, 0.8),
    layout = layout,
    topmargin = 10Plots.mm,
)
##

default(topmargin = 10Plots.mm)
get_data(syn) = load(datadir("up_down", "bistability", "balance_$syn.jld2"))
NMDA_balance = get_data("NMDA")["data"]
ratios = get_data("NMDA")["ratios"]
models = collect(values(NMDA_balance))
titles = ["distal\ndistal", "distal\nproximal", "proximal\nproximal", "soma\nonly"]

# models[1].critical
# heatmap(TN.νs, ratios, models[1].critical[2,..]', xscale=:log, c=:amp)
# hline!([1], c=:white)
# heatmap(TN.νs,ratios, NMDA_balance["distal-distal"].critical', xscale=:log, c=:amp)

models[1].critical[1, :, :]
plots_balance = []
my_greys = reverse(collect(palette(:greys, 5)))

for (title, (k, v)) in zip(titles, NMDA_balance)
    p = contour(
        TN.νs,
        ratios,
        reverse(v.critical[3, :, :]', dims = 1),
        clims = (0, 40),
        c = my_greys,
        cbar = false,
        ylabel = "",
        yticks = :none,
        xticks = :none,
        title = title,
        levels = 5,
        fill = true,
        lc = :greys,
    )
    push!(plots_balance, p)
    hline!([1], c = :black, ls = :dash, label = "", alpha = 0.2)
    # @show m
end

yticks = ([0.1, 1, 2], ["-", L"\bar k_{EI}", "+"])
plot!(plots_balance[1], yticks = yticks, yrotation = 0)
plot!(
    plots_balance[3],
    yticks = yticks,
    yrotation = 0,
    xticks = :true,
    ylabel = "                      E/I ratio",
    xlabel = "                                 Input rates (kHz)",
)
plot!(plots_balance[4], xticks = :true)
z3 = plot(plots_balance..., xscale = :log)

plots_membrane = []
for (title, (k, v)) in zip(titles, NMDA_balance)
    p = heatmap(
        TN.νs,
        ratios,
        reverse(v.μ[1, :, :]', dims = 1),
        c = :bluesreds,
        cbar = false,
        ylabel = "",
        yticks = :none,
        xticks = :none,
        title = title,
        clims = (-70, -50),
    )
    push!(plots_membrane, p)
    hline!([1], c = :black, ls = :dash, label = "", alpha = 0.2)
    # @show m
end

yticks = ([0.1, 1, 2], ["-", L"\bar k_{EI}", "+"])
plot!(plots_membrane[1], yticks = yticks, yrotation = 0)
plot!(
    plots_membrane[3],
    yticks = yticks,
    yrotation = 0,
    xticks = :true,
    ylabel = "                      E/I ratio",
    xlabel = "                                 Input rates (kHz)",
)
plot!(plots_membrane[4], xticks = :true)
z4 = plot(plots_membrane..., xscale = :log)

layout = @layout [a{0.40w} _ c{0.03w} _ a{0.40w}]
# path = mkpath(plotsdir("up_down","up_down","robustness"))
# savefig(z, joinpath(path,"critical_window.pdf"))
pp1 = colorbar(
    -65,
    -50,
    c = :bluesreds,
    label = "Membrane potential (mV)",
    ylabel = "Average membrane \n potential (mv) ",
    topmargin = 10Plots.mm,
    xlabel = " ",
)
pp2 = colorbar(
    0,
    50,
    c = :amp,
    label = "Membrane potential (mV)",
    ylabel = "Max bimodal \nwindow",
)
z5 = plot(pp2, pp1, layout = (2, 1))
zz2 = plot(
    z3,
    z5,
    z4,
    size = two_columns_size .* (1.2, 0.8),
    layout = layout,
    topmargin = 10Plots.mm,
    rightmargin = 5Plots.mm,
    bottommargin = 10mm,
    leftmargin = 5Plots.mm,
)
zz = plot(zz1, zz2, layout = (2, 1), size = two_columns_size .* (1.2, 1.5))

path = mkpath(plotsdir("up_down", "Figures"))
savefig(zz, joinpath(path, "Fig6ABC.pdf"))
