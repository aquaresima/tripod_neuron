using DrWatson
@quickactivate "Tripod"
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

EIs = string.(round.(collect([0.01, 0.25, 0.5, 0.8, 1]) .^ (-1), digits = 2))
EIs[1] = "∞"
CList = palette([:orange, :lightblue], length(ds))
plots = []
titles = ["NMDA+AMPA", "AMPA only", "", ""]
for (name, t) in zip([:dend_NMDA, :dend_AMPA, :istdp_NMDA, :istdp_AMPA], titles)
    p = plot(title = t)
    _data = getfield(opt_kies, name)
    for (l, c) in zip(ds[1:5:end], CList[5:5:end])
        d = findfirst(x -> x == l, ds)
        if d !== 1
            # && (c=:black)
            Plots.scatter!(
                p,
                νs,
                _data[:, d],
                xaxis = :log,
                markersize = 4,
                msc = c,
                c = c,
                xlims = (0.5, 55),
                label = "",
            )
        end
        @show ds[d]
    end
    # Plots.hline!([1], c=:black, ls=:dash)
    Plots.plot!(
        p,
        legend = false,
        xlabel = "                           Input rate (kHz)",
        ylabel = "E/I ratio",
    )
    Plots.push!(plots, p)
    Plots.scatter!(
        p,
        νs,
        _data[:, 1],
        xaxis = :log,
        markersize = 5,
        shape = :diamond,
        msc = :black,
        c = :black,
        xlims = (0.5, 80),
        label = "Soma",
        ylims = (-0.1, 1.2),
    )
    _yticks = ((0.01, 0.25, 0.5, 0.8, 1), EIs)
    plot!(yticks = _yticks)
end

_xticks_NMDA =
    ([νs[min_NMDA], 2, 5, 15, 50], [round.(νs[[min_NMDA]], digits = 1)..., 2.0, 5, 15, 50])
_xticks_AMPA =
    ([νs[min_AMPA], 2, 5, 15, 50], [round.(νs[[min_AMPA]], digits = 1)..., 2.0, 5, 15, 50])

plot!(plots[1], xlabel = "")
plot!(xticks = _xticks_NMDA, xrotation = 45)
Plots.vline!([νs[min_NMDA]], c = :black, ls = :dash)

plot!(plots[3])
Plots.vline!([νs[min_NMDA]], c = :black, ls = :dash)
plot!(xticks = _xticks_NMDA, xrotation = 45)

plot!(plots[2], xlabel = "", ylabel = "", yticks = :none, yaxis = false)
plot!(xticks = _xticks_AMPA, xrotation = 45)
plot!(twinx(), frame = :box, yticks = :none, ylabel = "Rate")
Plots.vline!([νs[min_AMPA]], c = :black, ls = :dash)

plot!(plots[4], ylabel = "", yticks = :none, xlabel = " ", yaxis = false)
plot!(xticks = _xticks_AMPA, xrotation = 45)
plot!(twinx(), frame = :box, yticks = :none, ylabel = "iSTDP")
Plots.vline!([νs[min_AMPA]], c = :black, ls = :dash)


hm = zeros(400, 1)
c = cgrad([:orange, :lightblue], length(ds))
hm[:, 1] .= collect(101:500)
# inset = (1, bbox(-0.9, 1.1, 1, 1, :bottom, :right))
xs = round.(Int, [100, 300, 500])
_xticks = ((0, 200, 400), xs)
plot(frame = :none, ticks = :none)
ss1 = heatmap!(
    twinx(),
    hm[:, :],
    c = c,
    cbar = false,
    yticks = _xticks,
    xticks = :none,
    frame = :box,
    ylabel = "Dendritic length (μm)",
    guidefontsize = 18,
    topmargin = 25Plots.mm,
    bottommargin = 25Plots.mm,
    rotation = -90,
)
# #
layout = @layout [grid(2, 1) grid(2, 1)]
p = plot(
    plots[1],
    plots[3],
    plots[2],
    plots[4],
    layout = layout,
    size = one_column_size,
    frame = :box,
)
path = mkpath(plotsdir("up_down", "balance"))
savefig(p, joinpath(path, "kie_values_50.pdf"))
plot(p, size = two_columns_size)
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
##

dd = ((data.dend_NMDA[:, 2:end] ./ mean(data.dend_NMDA[:, 2:end], dims = 2)) .- 1) .* 100
CList_b = collect(palette(:blues, 41, rev = true))
p1 = plot(
    νs,
    dd,
    c = CList',
    label = "NMDA",
    lw = 3,
    legend = false,
    xscale = :log,
    ylims = (-10, 10),
)
# p1= plot(νs, data.dend_NMDA[:,2:end], c=CList', label="NMDA", lw=3, legend=false, xscale=:log)
dd = (data.dend_AMPA[:, 2:end] ./ mean(data.dend_AMPA[:, 2:end], dims = 2) .- 1) .* 100
# dd = data.dend_AMPA[:,2:end] ./  mean(data.dend_AMPA[:,2:end], dims=2)-1
vline!([νs[min_NMDA]], c = :black, lw = 3, ls = :dash)
p2 = plot(
    νs,
    dd,
    c = CList',
    label = "AMP",
    lw = 3,
    legend = false,
    xscale = :log,
    ylims = (-100, 100),
)
dd = data.dend_AMPA[:, 1] ./ mean(data.dend_AMPA[:, 1], dims = 2)
vline!([νs[min_AMPA]], c = :black, lw = 3, ls = :dash)
# p3= plot(νs, dd,  c=:black, label="soma ", yticks=:none, ylims=(0.1,5), yscale=:log)
plot!(p1, yaxis = "E/I ratio (rate)")
plot!(
    p2,
    xlabel = "Input rate (kHz)",
    legend = :topleft,
    bottommargin = 10Plots.mm,
    yticks = ([-100, 0, 100], ["-100%", "0%", "+100%"]),
    title = "AMPA only",
)
plot!(
    p1,
    xlabel = "Input rate (kHz)",
    legend = :topleft,
    bottommargin = 10Plots.mm,
    yticks = ([-10, 0, 10], ["-10%", "0%", "+10%"]),
    title = "NMDA+AMPA",
)
q2 = plot(
    p1,
    p2,
    legend = false,
    layout = (1, 2),
    xscale = :log10,
    leftmargin = 0Plots.mm,
    rightmargin = 0Plots.mm,
)
##
# q2= plot!(νs, mean(data.dend_NMDA[:,2:end], dims=2), ribbon=std(data.dend_NMDA, dims=2), c=:darkblue, label="NMDA", lw=3)
# q2= plot!(νs, mean(data.dend_AMPA[:,2:end], dims=2), ribbon=std(data.dend_AMPA, dims=2), c=:darkred, label="AMPA", lw=3, yaxis=L"\langle"*"E/I ratio"* L"\rangle", xscale=:log10, xlabel="Input rate (kHz)", legend=:topleft, bottommargin=10Plots.mm)
# plot!(yticks=_yticks, ylims=(0,1))
# plot!(xticks=([νs[13], 1.5,5,15,50], [round(νs[13], digits=1), 1.5,5,15,50]))


q = plot(q1, q2, layout = (2, 1), leftmargin = 5Plots.mm)
plot(q, size = (800, 700), leftmargin = 5Plots.mm)
##
# savefig(p, joinpath(path,"kie_cor.pdf"))


get_data(title, type = "data") =
    load(datadir("up_down", "activity", "full_activity_$(title)_highres.jld2"))
AMPA = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)
νs = collect(exp.(range(log(0.1), log(5), length = 30)))
@unpack min_AMPA, min_NMDA = TN

sum(AMPA.mean_v[1, :, :, :])
c = collect(palette(:roma, size(AMPA.mean_v)[end]...))'
layout = @layout[[aa{0.01w} a{0.460w} b{0.46w}]; [aa{0.01w} a{0.46w} b{0.46w}]]
# plot!(layout=layout)


default(legend = false, margin = 5Plots.mm)
NMDA.mean_v
my_diag(NMDA.mean_v, 2)
NMDA.mean_v
_xticks = ([100, 300, 500], [100, 300, 500])
_yticks = (-70:15:-20, -70:15:-20)

default(bottommargin = 0Plots.mm, rightmargin = 0Plots.mm)
s1 = plot(
    ds[2:end],
    my_diag(NMDA.mean_v, 2)[2:end, :],
    title = "AMPA + NMDA",
    c = c,
    cbar = false,
    xticks = :none,
    yticks = _yticks,
    legend = false,
    ylims = (-70, -20),
)
plot!(ds[2:end], my_diag(NMDA.mean_v, 2)[2:end, min_NMDA], title = "", c = :white)
plot!(
    ds[2:end],
    my_diag(NMDA.mean_v, 2)[2:end, min_NMDA],
    title = "",
    c = :black,
    ls = :dash,
)
s2 = plot(
    ds[2:end],
    my_diag(AMPA.mean_v, 2)[2:end, :],
    title = "AMPA only",
    c = c,
    cbar = false,
    xticks = :none,
    yticks = :none,
    legend = false,
    ylims = (-70, -20),
)
plot!(twinx(), frame = :box, yticks = :none, ylabel = "Dendrite")
s4 = plot!(ds[2:end], my_diag(AMPA.mean_v, 2)[2:end, min_AMPA], title = "", c = :white)
s4 = plot!(
    ds[2:end],
    my_diag(AMPA.mean_v, 2)[2:end, min_AMPA],
    title = "",
    c = :black,
    ls = :dash,
)


default(bottommargin = 5Plots.mm, rightmargin = 0Plots.mm)
s3 = plot(
    ds[2:end],
    my_diag(NMDA.mean_v)[2:end, :],
    title = "",
    c = c,
    cbar = false,
    ylabel = "                     Membrane potential (mV)",
    xlabel = "                          Dendritic length (μm)",
    xticks = _xticks,
    yticks = (-65:10:-55, -65:10:-55),
    ylims = (-70, -50),
    legend = false,
)
plot!(ds[2:end], my_diag(NMDA.mean_v)[2:end, min_NMDA], title = "", c = :white)
plot!(ds[2:end], my_diag(NMDA.mean_v)[2:end, min_NMDA], title = "", c = :black, ls = :dash)
s4 = plot(
    ds[2:end],
    my_diag(AMPA.mean_v)[2:end, :],
    title = "",
    c = c,
    cbar = false,
    yticks = :none,
    xlabel = " ",
    xticks = _xticks,
    guidefontsize = 18,
    ylims = (-70, -50),
    legend = false,
)
plot!(twinx(), yticks = :none, ylabel = "Soma")

s4 = plot!(ds[2:end], my_diag(AMPA.mean_v)[2:end, min_AMPA], title = "", c = :white)
s4 = plot!(
    ds[2:end],
    my_diag(AMPA.mean_v)[2:end, min_AMPA],
    title = "",
    c = :black,
    ls = :dash,
)


hm = zeros(100, 1)
hm[:, 1] .= collect(1:length(hm))
inset = (1, bbox(-0.9, 1.1, 1, 1, :bottom, :right))
xs = round.(Int, range(1, 100, 3))
xxs = round.(νs[[10, 20, 30]], digits = 1)
_xticks = (xs, xxs)
plot(frame = :none)
ss2 = heatmap!(
    twinx(),
    hm[1:100, :],
    c = :roma,
    cbar = false,
    xticks = :none,
    frame = :box,
    ylabel = "Input rate (kHz)",
    guidefontsize = 15,
    yticks = _xticks,
    topmargin = 15Plots.mm,
    bottommargin = 15Plots.mm,
    rotation = -90,
    rightmargin = 10Plots.mm,
)

layout = @layout [grid(2, 1) grid(2, 1)]# b{0.05w}] 
s = plot(
    s1,
    s3,
    s2,
    s4,
    legend = false,
    layout = layout,
    frame = :box,
    size = one_column_size .* (1.0, 1.5),
    xrotation = 45,
)
# layout = @layout [a{0.95w} b{0.05w}]




##
layout = @layout [
    a{0.45w} _ a{0.45w} z{0.03w} _
    # _ _ _
    a{0.45w} z{0.03w} _ a{0.45w}
]
z = plot(
    plot(frame = :none),
    p,
    ss1,
    s,
    ss2,
    q,
    size = two_columns_size .* (1.25, 1.7),
    layout = layout,
    rightmargin = 0Plots.mm,
)
path = mkpath(plotsdir("up_down", "Figures"))
savefig(z, joinpath(path, "Fig3.pdf"))
z
##


## Further analysis
CList = palette([:orange, :lightblue], length(ds))
plot()
for x = 1:42
    plot!(data.dend_NMDA[:, x], data.istdp_NMDA[:, x], c = CList[x], label = ds[x], lw = 2)
end
plot!()


plot(
    mean(NMDA.rate[1, 2:end, 2:end, :], dims = [1, 2])[1, 1, :],
    rightmargin = 20mm,
    colorbar = true,
)
plot!(
    mean(AMPA.rate[1, 2:end, 2:end, :], dims = [1, 2])[1, 1, :],
    rightmargin = 20mm,
    colorbar = true,
)
plot!(AMPA.rate[1, 1, 1, :])
hline!([-55])

##
d1, d2 = 20, 20
nmda = true
model = (ds = (TN.ds[d1], TN.ds[d2]), nmda = nmda, soma_only = sum(TN.ds[d1]) == 0)
rate = zeros(44)
simtime = 50_000
for f = 1:44
    cond = TN.get_balance_conditions(d1, d2, f, nmda = nmda)
    inputs = TN.make_spikes(simtime; cond...)
    v_spikes = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
    rate[f] = TN.get_spike_rate(v_spikes)
end

plot(rate)


simtime = 100000
model = TN.models[1]
cond = TN.get_balance_conditions(model.ds..., 20, nmda = true)
β = 600
inputs = TN.make_spikes(simtime, β; cond...)
voltage = TN.run_tripod(
    inputs,
    simtime;
    model...,
    AdEx = TN.AdEx,
    syn_model = cond.syn_model,
    do_spikes = true,
)
plot(voltage[1, end-10000:end])

histogram(
    voltage[1, :],
    bins = 1000,
    normed = true,
    xlims = (-70, -50),
    ylims = (0, 0.3),
    xlabel = "Membrane potential (mV)",
    ylabel = "Probability density",
    legend = false,
    frame = :box,
    size = (600, 500),
)
