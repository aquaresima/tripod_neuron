using DrWatson
@quickactivate "Tripod"
using Plots, JLD2, HDF5, Statistics
using Plots, JLD2, HDF5, Statistics
using Revise
using TripodNeuron
include(projectdir("scripts", "plots", "default_plots.jl"))
## Load data from previous computations

get_data(title, type = "data"; suffix = "") =
    load(datadir("up_down", "balance", "inhibitory_balance_$title$suffix.jld2"))[type]

μmem = -55.0
data = load(datadir("up_down", "balance", "optimal_IE_μmem$(μmem).jld2"))

@unpack min_AMPA, min_NMDA, νs, opt_kies = data
@unpack ds = TN
##
data = opt_kies
dend_AMPA = data.dend_AMPA[min_AMPA:end, 2:end]
dend_NMDA = data.dend_NMDA[min_NMDA:end, 2:end]
istdp_AMPA = data.istdp_AMPA[min_AMPA:end, 2:end]
istdp_NMDA = data.istdp_NMDA[min_NMDA:end, 2:end]
istdp_soma = (data.istdp_AMPA[min_AMPA:end, 1] .+ data.istdp_NMDA[min_AMPA:end, 1]) ./ 2
dend_soma = (data.dend_AMPA[min_AMPA:end, 1] .+ data.dend_NMDA[min_AMPA:end, 1]) ./ 2
q1 = plot(topmargin = 13mm)

##
using StatsBase
using LsqFit
using LaTeXStrings

data = opt_kies
dend_AMPA = data.dend_AMPA[min_AMPA:end, 2:end]
dend_NMDA = data.dend_NMDA[min_NMDA:end, 2:end]
istdp_AMPA = data.istdp_AMPA[min_AMPA:end, 2:end]
istdp_NMDA = data.istdp_NMDA[min_NMDA:end, 2:end]
istdp_soma = (data.istdp_AMPA[min_AMPA:end, 1] .+ data.istdp_NMDA[min_AMPA:end, 1]) ./ 2
dend_soma = (data.dend_AMPA[min_AMPA:end, 1] .+ data.dend_NMDA[min_AMPA:end, 1]) ./ 2
q1 = plot(topmargin = 13mm)
CList = collect(palette([:orange, :lightblue], length(ds)))
# CList_r = collect(palette(:reds, 41, rev=true))

p1 = plot(yticks = :none)
for r in collect(1:length(νs[min_AMPA:end]))[1:3:39]
    Plots.scatter!(istdp_AMPA[r, :], dend_AMPA[r, :], c = CList, msc = CList)
    plot!([0, 1], [0, 1], c = :black, lw = 3, ls = :dash)
end
plot(p1, legend = false)
p2 = plot(yticks = :none)
for r in collect(1:length(νs[min_NMDA:end]))[1:3:39]
    Plots.scatter!(istdp_NMDA[r, :], dend_NMDA[r, :], c = CList, msc = CList)
    plot!([0, 1], [0, 1], c = :black, lw = 3, ls = :dash)
end
p3 = plot(yticks = :none)
for r = 1:length(νs[min_AMPA:end])
    Plots.scatter!(istdp_soma[r, :], dend_soma[r, :], c = :black, msc = :black)
    plot!([0, 1], [0, 1], c = :black, lw = 3, ls = :dash)
end

_yticks = ((0.01, 0.25, 0.5, 0.8, 1), EIs)
plot!(
    p2,
    ylabel = "E/I ratio (rate)",
    alpha = 1.0,
    legend = false,
    yticks = ([0, 0.5, 1], [0, 0.5, 1]),
    ylims = (0, 1),
    xlims = (0, 1),
)
plot!(yticks = _yticks)
plot!(p1, xlabel = "E/I ratio (iSTDP)")
EIs = string.(round.(collect([0.01, 0.25, 0.5, 0.8, 1]) .^ (-1), digits = 2))
EIs[1] = "∞"
# plot!(p2,xticks=_yticks)
# plot!(xticks=_yticks)
@. linear_model(x, p) = x * p[1] + p[2]

nmda_fit = curve_fit(linear_model, istdp_NMDA[:], dend_NMDA[:], [0.0, 0.0])
xx = [minimum(istdp_NMDA), maximum(istdp_NMDA)]
yy = xx .* nmda_fit.param[1] .+ nmda_fit.param[2]
ampa_fit = curve_fit(linear_model, istdp_AMPA[:], dend_AMPA[:], [0.0, 0.0])
xx = [minimum(istdp_AMPA), maximum(istdp_AMPA)]
yy = xx .* ampa_fit.param[1] .+ ampa_fit.param[2]
c = round(cor(istdp_NMDA[:], dend_NMDA[:]), digits = 2)
Plots.annotate!(
    p2,
    [(0.3, 0.10, Plots.text(L"\rho" * " = $c", :black, 14, rotation = 0, :left))],
)
plot!(p2, title = "NMDA+AMPA")
c = round(cor(istdp_AMPA[:], dend_AMPA[:]), digits = 2)
Plots.annotate!(
    p1,
    [(0.3, 0.1, Plots.text(L"\rho" * " = $c", :black, 14, rotation = 0, :left))],
)
plot!(p1, title = "AMPA only")
c = round(cor(istdp_soma[:], dend_soma[:]), digits = 2)
Plots.annotate!(
    p3,
    [(0.3, 0.1, Plots.text(L"\rho  " * " = $c", :black, 14, rotation = 0, :left))],
)
plot!(p3, title = "Soma only")
q1 = plot(
    p2,
    p1,
    p3,
    layout = (1, 3),
    xrotation = 90,
    leftmargin = 0Plots.mm,
    rightmargin = 0Plots.mm,
    legend = false,
)
plot!(xticks = _yticks)
p = plot!(
    size = (800, 400),
    leftmargin = 10Plots.mm,
    legend = false,
    bottommargin = 15Plots.mm,
)
savefig(p, plotsdir("up_down", "Figures", "Fig10_Appendix.pdf"))
