##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
using EllipsisNotation

include(projectdir("scripts", "plots", "default_plots.jl"))

##
get_data(title, type = "data") =
    load(datadir("up_down", "asymmetries", "rate_β=0$title.jld2"))
AMPA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)

@unpack νs, ds = TN
v = NMDA_as.cors

default(lw = 4)
p = plot()
plot!(
    νs,
    mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    label = "All",
    ls = :dash,
    c = :black,
)
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
vline!([νs[TN.min_NMDA]], c = :black, ls = :dash, label = "AMPA")
plot!(legend = :false, ylims = (0, 1), title = "NMDA+AMPA")

v = AMPA_as.cors
q = plot()
plot!(
    νs,
    mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    label = "All",
    ls = :dash,
    c = :black,
)
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
vline!([νs[TN.min_AMPA]], c = :black, ls = :dash, label = "AMPA")
plot!(legend = :topright, ylims = (0, 1), title = "AMPA only")
#
p1 = plot(
    p,
    q,
    layout = (1, 2),
    size = (800, 400),
    ylabel = "Dendritic correlation",
    xlabel = "Input rate (kHz)",
    xscale = :log10,
)
## Rate
v = NMDA.rate
TN.min_NMDA

p = plot()
# plot!(νs, mean(v[1, 1:1, 1:1, :], dims = (1, 2))[1, 1, :], label = "Soma")
plot!(νs, mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All")
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :false, ylims = (0, 5), title = "NMDA + AMPA")
annotate!(0.15, 6.2, text("A", 20, :center))

v = AMPA.rate
q = plot()
# plot!(νs, mean(v[1, 1:1, 1:1, :], dims = (1, 2))[1, 1, :], label = "Soma")
plot!(νs, mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All")
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :topright, ylims = (0, 5), title = "AMPA only")
#
p2 = plot(
    p,
    q,
    layout = (1, 2),
    size = (800, 400),
    ylabel = "Firing rate (Hz)",
    xlabel = "Input rate (kHz)",
    xscale = :log10,
)
## ISI CV
v = NMDA.rate

p = plot()
plot!(νs, mean(v[2, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All", c=:black, ls=:dash)
plot!(νs, mean(v[2, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[2, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[2, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :false, ylims = (0, 2), title = "NMDA + AMPA")
annotate!(0.1, 2.2, text("B", 20, :center))

v = AMPA.rate
q = plot()
plot!(νs, mean(v[2, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All", c=:black, ls=:dash)
plot!(νs, mean(v[2, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[2, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[2, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :false, ylims = (0, 2), title = "AMPA only")
#
p3 = plot(
    p,
    q,
    layout = (1, 2),
    size = (800, 400),
    ylabel = "ISI CV",
    xlabel = "Input rate (kHz)",
    xscale = :log10,
)

## STD
get_data(title, type = "data") =
    load(datadir("up_down", "activity", "full_activity_$title.jld2"))
AMPA = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)
v = NMDA.std_v

p = plot()
plot!(νs, mean(v[1, 1:1, 1:1, :], dims = (1, 2))[1, 1, :], label = "Soma")
plot!(νs, mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All")
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :false, ylims = (0, 5), title = "NMDA + AMPA")
annotate!(0.15, 6.2, text("C", 20, :center))

v = AMPA.std_v
q = plot()
plot!(νs, mean(v[1, 1:1, 1:1, :], dims = (1, 2))[1, 1, :], label = "Soma")
plot!(νs, mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :], label = "All")
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(legend = :false, ylims = (0, 5), title = "AMPA only")
p4 = plot(
    p,
    q,
    layout = (1, 2),
    size = (800, 400),
    ylabel = "Mem. std (mV) ",
    xlabel = "Input rate (kHz)",
    xscale = :log10,
)
#
p = plot(p2, p3, p4, p1, layout = (4, 1), size = (800, 1200), leftmargin = 10mm)
p = plot(p4,p3, p1, layout=(3,1), size=(800,1000))

savefig(p, plotsdir("up_down", "Figures", "Fig7_Appendix.pdf"))
#:])



get_data(title, type = "data") =
    load(datadir("up_down", "asymmetries", "rate_β=0$title.jld2"))
AMPA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)

@unpack νs, ds = TN

v = NMDA_as.cors

default(lw = 4)
qq2 = plot()
# plot!(νs, mean(v[1,1:1,1:1,:], dims=(1,2))[1,1,:], label="Soma")
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(
    νs,
    mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    label = "All",
    ls = :dot,
    c = :black,
)
vline!([νs[TN.min_NMDA]], c = :black, ls = :dash, lw = 3)
plot!(
    legend = :false,
    ylims = (0, 1),
    title = "NMDA+AMPA",
    ylabel = "Dendritic correlation",
)

v = AMPA_as.cors
qq1 = plot()
plot!(νs, mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :], label = "Proximal")
plot!(νs, mean(v[1, 23:end, 23:end, :], dims = (1, 2))[1, 1, :], label = "Distal")
plot!(νs, mean(v[1, 23:end, 2:22, :], dims = (1, 2))[1, 1, :], label = "Asymmetrical")
plot!(
    νs,
    mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    label = "All",
    ls = :dot,
    c = :black,
)
vline!([νs[TN.min_AMPA]], c = :black, ls = :dash, label = "", lw = 3)
plot!(legend = :topright, ylims = (0, 1), title = "AMPA only")
#
q1 = plot(
    qq2,
    qq1,
    layout = (1, 2),
    size = (800, 400),
    xlabel = "Input rate (kHz)",
    xscale = :log10,
    leftmargin = 10Plots.mm,
)

# plot(q, size=(800, 700))
#

get_data(title, type = "data") =
    load(datadir("up_down", "asymmetries", "rate_β=0$title.jld2"))
AMPA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA_as = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)

@unpack νs, ds = TN
# ampa_color= :darkred
# nmda_color= :darkblue
# ## correlation

v = NMDA_as.cors

default(lw = 4)
qq2 = plot()
# plot!(νs, mean(v[1,1:1,1:1,:], dims=(1,2))[1,1,:], label="Soma")
plot!(νs, mean(v[2:3, 2:22, 2:22, :], dims = (1, 2, 3))[1, 1, 1, :], label = "Proximal")
plot!(νs, mean(v[2:3, 23:end, 23:end, :], dims = (1, 2, 3))[1, 1, 1, :], label = "Distal")
plot!(
    νs,
    mean(v[2:3, 23:end, 2:22, :], dims = (1, 2, 3))[1, 1, 1, :],
    label = "Asymmetrical",
)
plot!(
    νs,
    mean(v[2:3, 2:end, 2:end, :], dims = (1, 2, 3))[1, 1, 1, :],
    label = "All",
    ls = :dot,
    c = :black,
)
vline!([νs[TN.min_NMDA]], c = :black, ls = :dash, lw = 3)
plot!(
    legend = :false,
    ylims = (0, 1),
    title = "NMDA+AMPA",
    ylabel = "Soma - dendrites\ncorrelation",
)

v = AMPA_as.cors
qq1 = plot()
plot!(νs, mean(v[2:3, 2:22, 2:22, :], dims = (1, 2, 3))[1, 1, 1, :], label = "Proximal")
plot!(νs, mean(v[2:3, 23:end, 23:end, :], dims = (1, 2, 3))[1, 1, 1, :], label = "Distal")
plot!(
    νs,
    mean(v[2:3, 23:end, 2:22, :], dims = (1, 2, 3))[1, 1, 1, :],
    label = "Asymmetrical",
)
plot!(
    νs,
    mean(v[2:3, 2:end, 2:end, :], dims = (1, 2, 3))[1, 1, 1, :],
    label = "All",
    ls = :dot,
    c = :black,
)
vline!([νs[TN.min_AMPA]], c = :black, ls = :dash, label = "", lw = 3)
plot!(legend = :topright, ylims = (0, 1), title = "AMPA only")
#
q2 = plot(
    qq2,
    qq1,
    layout = (1, 2),
    size = (800, 400),
    xlabel = "Input rate (kHz)",
    xscale = :log10,
    leftmargin = 10Plots.mm,
)

p = plot(
    q1,
    q2,
    size = (1300, 320),
    layout = (1, 2),
    bottommargin = 20Plots.mm,
    leftmargin = 20Plots.mm,
)
savefig(p, plotsdir("up_down", "Figures", "Fig11_Appendix.pdf"))
# plot(q, size=(800, 700))z
p
