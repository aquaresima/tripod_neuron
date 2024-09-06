using Pkg
Pkg.activate("@v1.10")
using StatsPlots, UncertainData, RollingFunctions

using DrWatson
@quickactivate "Tripod"
using TripodNeuron

using Plots, Revise, Random
using LsqFit
using Statistics
using HypothesisTests
using EffectSizes

include(projectdir("scripts", "plots", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_analysis.jl"))
##

path = mkpath(plotsdir("up-down", "Figures"))
simtime = 1_000_000
β = 200
τ = 50
interval = 150_000:200_000
@unpack νs = TN
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
_rate = 11;
@show νs[_rate];
##
model = TN.models[3]
Random.seed!(133)
cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
voltage_NMDA =
    TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
plot(voltage_NMDA[1, interval])
plot(cap_voltage(voltage_NMDA[1, interval]))
# histogram(voltage[1,:], bins=-80:-50, )
NMDA = updown_analysis(voltage_NMDA[1, :])
Random.seed!(133)
cond = TN.get_balance_conditions(model.ds..., _rate; nmda = false)
inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
voltage_AMPA =
    TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
# histogram(voltage[1,:], bins=-80:-50, )
AMPA = updown_analysis(voltage_AMPA[1, :])

model = TN.models[4]
Random.seed!(133)
cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
voltage_soma =
    TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
# histogram(voltage[1,:], bins=-80:-50, )
SOMA = updown_analysis(voltage_soma[1, :])


##


p1 = plot(
    plot(voltage_NMDA[1, 1:100_001], xticks = :none, c = :darkblue, xlabel = ""),
    plot(
        voltage_AMPA[1, 1:100_001],
        xticks = :none,
        c = :darkred,
        xlabel = "",
        ylabel = "Mem. potential (mV)",
    ),
    plot(
        voltage_soma[1, 1:100_001],
        c = :black,
        xlabel = "Time (s)",
        xticks = (0:20000:100_000, 0:2:10),
        bottommargin = 15mm,
    ),
    layout = (3, 1),
    legend = false,
    ylims = (-70, -30),
    margin = 5Plots.mm,
    yticks = (-70:20:-30, ["-70", "-50", "-30"]),
)
##
p1

# states = NMDA.up[1,:] ./ 1000
# _hist = length.(UncertainData.bin(0:0.01:10,states,1:length(states)))
# plot(0.01:0.01:10,cum_prob(_hist))
# states = NMDA.down[1,:] ./ 1000
# _hist = length.(UncertainData.bin(0:0.01:10,states,1:length(states)))
# plot!(0.01:0.01:10,cum_prob(_hist))
## Spike voltage correlation
f(x) = x > 0.001 ? "p = $(round(x, digits=2))" : "p < 10ˆ-3 "
nonzero = findall(NMDA.up[4, :] .> 0)
cNMDA = cor(NMDA.up[2, nonzero], NMDA.up[4, nonzero]) |> x->round(x, digits=2)
pnmda = pvalue(OneSampleZTest(atanh(cNMDA),1, length(NMDA.up[2,nonzero]))) |> f
nonzero = findall(AMPA.up[4, :] .> 0)
cAMPA = cor(AMPA.up[2, nonzero], AMPA.up[4, nonzero]) |> x->round(x, digits=2)
pampa = pvalue(OneSampleZTest(atanh(cAMPA),1, length(AMPA.up[2,nonzero]))) |> f
nonzero = findall(SOMA.up[4, :] .> 0)
cSOMA = cor(SOMA.up[2, nonzero], SOMA.up[4, nonzero])
psoma = pvalue(OneSampleZTest(atanh(cSOMA),1, length(SOMA.up[2,nonzero]))) |>f
nonzero = findall(SOMA.up[4, :] .> 0)

cortest(x,y) =
    if length(x) == length(y)
        2 * ccdf(Normal(), atanh(abs(cor(x, y))) * sqrt(length(x) - 3))
    else
        error("x and y have different lengths")
    end

using HypothesisTests

NMDA.up[2,nonzero]



p2 = plot(
    scatter(
        NMDA.up[2, :],
        NMDA.up[4, :],
        ylabel = "Firing rate (Hz)",
        c = :darkblue,
        title = "NMDA \nρ = $cNMDA; $pnmda",
        msc = :darkblue,
        smooth=true,
        lc=:grey,
        lw = 4,
        ms = 4,
    ),
    scatter(
        AMPA.up[2, :],
        AMPA.up[4, :],
        xlabel = "Avg. membrane potential (mV)",
        c = :darkred,
        title = "AMPA \nρ = $cAMPA; $pampa",
        msc = :darkred,
        smooth=true,
        lc=:grey,
        lw = 4,
        ms = 4,
    ),
    # scatter(
    #     SOMA.up[2, :],
    #     SOMA.up[4, :],
    #     # line=true,
    #     c = :black,
    #     title = "SOMA, \np = $(round(psoma, digits=2))",
    #     msc = :black,
    #     lc=:red,
    #     ms = 4,
    # ),
    layout = (1, 2),
    xlims = (-64, -52),
    legend = false,
    # xrotation = -45,
    xticks = (-62:5:-50),
    yticks = 0:30:60,
    ylims = (-5, 60),
    titlefontsize=10,
    frame=:axes, 
    rightmargin=3Plots.mm,
    leftmargin=3Plots.mm
)
p2_intermedia = plot(p2,bottommargin = 10Plots.mm, size=(600,400))
# savefig(p2, joinpath(path, "Fig7_firing rate.pdf"))
p2
##

_std_err =
    std.([
        NMDA.up[3, :],
        NMDA.down[3, :],
        AMPA.up[3, :],
        AMPA.down[3, :],
        # SOMA.up[3, :],
        # SOMA.down[3, :],
    ])
_std =
    mean.([
        NMDA.up[3, :],
        NMDA.down[3, :],
        AMPA.up[3, :],
        AMPA.down[3, :],
        # SOMA.up[3, :],
        # SOMA.down[3, :],
    ])
_mean = ([
    NMDA.up[2, :],
    NMDA.down[2, :],
    AMPA.up[2, :],
    AMPA.down[2, :],
    # SOMA.up[2, :],
    # SOMA.down[2, :],
])
p3 = boxplot(
    _mean,
    legend = false,
    c = [:darkblue :darkblue :darkred :darkred :black :black],
    xticks = ([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN"]),
    xrotation = -45,
    lc = :match,
)
plot!(p3, ylabel = "Soma potential (mV)", xlabel = "State")
p4 = bar(
    1:2,
    _std[1:2],
    yerror = _std_err[1:2],
    legend = false,
    c = :darkblue,
    xticks = ([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN"]),
    xrotation = -45,
    lc = :match,
)
bar!(
    3:4,
    _std[3:4],
    yerror = _std_err[3:4],
    legend = false,
    c = :darkred,
    xticks = ([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN"]),
    xrotation = -45,
    lc = :match,
)
# bar!(
#     5:6,
#     _std[5:6],
#     yerror = _std_err[5:6],
#     legend = false,
#     c = :black,
#     xticks = ([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN", "UP", "DOWN"]),
#     xrotation = -45,
#     lc = :match,
# )
plot!(p4, ylabel = "STD Membrane\npotential (mV)", xlabel = "State")
plot!(bottommargin = 10Plots.mm)

## Check P values

p_test = Vector{Float64}(undef, 3)
for n in [1, 2]
    i = 2 * n - 1
    m = mean.(_mean)
    v = std.(_mean)
    l = round.(Int, length.(_mean) ) 
    # @info OneSampleTTest(m[i]-m[i+1], v[i], l[i])
    @info string(mean(_mean[i]), " ", mean(_mean[i+1]))
    cd = CohenD(_mean[i],_mean[i+1], quantile=0.99)
    ci = confint(cd)
    @info string(effectsize(cd),"-", lower(ci), "-", upper(ci))
end

_mean

std(_mean[1])
OneSampleTTest

# length.(_mean[5])
p_test
_std_err

##
layout = @layout [
    a{0.7h}
    b{0.3h}
]
pp1 =
    plot(p1, p2, layout = layout, size = (700, 800), legend = false, leftmargin = 5Plots.mm)
plot!(p3, ylims=(-70,-50))
pp2 = plot(
    p3,
    p4,
    layout = (1,2),
    size = (600, 400),
    legend = false,
    topmargin = 10mm,
    bottommarging = 10mm,
    frame= :axes
)

savefig(pp2, joinpath(path, "Fig7_experimental_average.pdf"))
pp = plot(
    pp1,
    pp2,
    layout = (1, 2),
    size = (1400, 800),
    legend = false,
    leftmargin = 15Plots.mm,
)
##

# data = load(datadir("up_down","bistability","characterize_updown_(150, 150).jld2")) |> dict2ntuple
data = load(datadir("up_down", "bistability", "stimulus_condition.jld2")) |> dict2ntuple
function plot_cum_prob(data, bs, state, receptor = "NMDA"; color = :black, kwargs...)
    plot()
    # c = palette(, length(bs))
    for (b, ls) in zip(bs, [:solid :dash :dot])
        states = getfield(data[b][receptor], state)[1, :] ./ 1000
        _hist = length.(UncertainData.bin(0:0.01:10, states, 1:length(states)))
        plot!(0.01:0.01:10, cum_prob(_hist), c = color, ls = ls, label = "")
    end
    plot!(; kwargs...)
end

# βs[5]
default(grid = false)
q = plot_cum_prob(
    data.condition,
    [1, 6],
    :up,
    "SOMA",
    color = :black,
    legend = false,
    yticks = :none,
    rightmargin = 10Plots.mm,
)
q = plot!(twinx(), ylabel = "UP", title = "SOMA", yticks = :none, frame = :box, grid = true)
plot!(
    [[], []],
    ls = [:solid, :dash],
    label = ["ISI CV = 1." "ISI CV = 1.2"],
    c = :black,
    lw = 2,
)

z1 = plot(
    plot_cum_prob(
        data.condition,
        [1, 6],
        :up,
        "NMDA",
        color = :darkblue,
        title = "NMDA",
        legend = false,
        leftmargin = 10Plots.mm,
        # yticks = (0:0.5:1, 0:0.5:1),
    ),
    # plot_cum_prob(
    #     data.condition,
    #     [1, 6],
    #     :up,
    #     "AMPA",
    #     color = :darkred,
    #     title = "AMPA",
    #     legend = false,
    #     yticks = :none,
    # ),
    # q,
    xticks = (1:2:8),
    ylims = (0, 1),
    xlims = (0, 8),
    topmargin = 10Plots.mm,
    yticks = :none,
    title= "UP"
    # xticks = :none,
)

q = plot_cum_prob(
    data.condition,
    [1, 6],
    :down,
    "SOMA",
    color = :black,
    legend = false,
    yticks = :none,
    rightmargin = 10Plots.mm,
)
q = plot!(twinx(), ylabel = "DOWN", yticks = :none, frame = :box)
z2 = plot(
    plot_cum_prob(
        data.condition,
        [1, 6],
        :down,
        "NMDA",
        color = :darkblue,
        legend = false,
        ylabel="Cumulative probability (%)",
        leftmargin = 10Plots.mm,
        yticks = (0:0.2:1, 0:20:100),
    ),
    # plot_cum_prob(
    #     data.condition,
    #     [1, 6],
    #     :down,
    #     "AMPA",
    #     color = :darkred,
    #     legend = false,
    #     yticks = :none,
    # ),
    # q,
    layout = (1, 3),
    legend = false,
    xlabel = "Time (s)",
    xticks = (1:2:5),
    ylims = (0, 1),
    xlims = (0, 5),
    bottommargin = 10Plots.mm,
    title= "DOWN",
    topmargin = 3mm,
)
zz1 = plot(z2, z1, layout = (1, 3), size = (400, 400), leftmargin=3Plots.mm, rightmargin=3Plots.mm, frame=:axes)
plot!(zz1, subplot=2, rightmargin=10Plots.mm)
plot!(zz1, subplot=1, leftmargin=10Plots.mm, size=(600,400))
##

## APPENDIX

plots = []
for (v, c, t) in zip(
    [
        cap_voltage(voltage_NMDA[1, :]),
        cap_voltage(voltage_AMPA[1, :]),
        # cap_voltage(voltage_soma[1, :]),
    ],
    [:darkblue, :darkred, :black],
    ["NMDA", "AMPA", "SOMA"],
)
    up, down, thresh_up, thresh_down = voltage_cross_counter(v, 3 / 4, 1 / 4)
    upstart = [(x[1]-1000:1:x[1]+1000) for x in up[3:end]]
    p = plot()
    [plot!(v[x], alpha = 0.05, c = :black) for x in upstart[1:10:end]]
    plot!(
        mean(hcat([v[x] for x in upstart]...), dims = 2)[:, 1],
        legend = false,
        c = c,
        lw=3,
        title = t,
    )
    plot!(yticks = :none)
    push!(plots, p)
end
plot!(
    plots[1],
    yticks = ([-75, -65, -55]),
    ylabel = " ",
    title = "NMDA",
    legend = false,
    leftmargin = 25Plots.mm,
)
qq1 = plot(
    plots...,
    layout = (1, 2),
    size = (700, 200),
    legend = false,
    leftmargin = 5Plots.mm,
    xticks = :none,
    ylims = (-75, -50),
    xlims = (0, 1600),
    topmargin = 10Plots.mm,
)
plots = []


v = cap_voltage(voltage_NMDA[1, :])
up, down, thresh_up, thresh_down = voltage_cross_counter(v, 3 / 4, 1 / 4)
up[3:end]
for (v, c, t) in zip(
    [
        cap_voltage(voltage_NMDA[1, :]),
        cap_voltage(voltage_AMPA[1, :]),
        # cap_voltage(voltage_soma[1, :]),
    ],
    [:darkblue, :darkred, :black],
    ["NMDA", "AMPA", "SOMA"],
)
    up, down, thresh_up, thresh_down = voltage_cross_counter(v, 3 / 4, 1 / 4)
    downstart = [(x[1]-1000:1:x[1]+1000) for x in down[3:end]]
    p = plot()
    [plot!(v[x], alpha = 0.05, c = :black) for x in downstart[1:10:end]]
    plot!(mean(hcat([v[x] for x in downstart]...), dims = 2)[:, 1], legend = false, c = c)
    plot!(yticks = :none)
    push!(plots, p)
end
plot!(
    plots[1],
    yticks = ([-75, -65, -55]),
    ylabel = "Membrane potential (mV)",
    legend = false,
    leftmargin = 15Plots.mm,
)
plot!(plots[2], xlabel = "Time to state change (ms)")
qq2 = plot(
    plots...,
    layout = (1, 2),
    size = (700, 200),
    legend = false,
    xticks = ([0, 1000, 2000], [-100, 0, 100]),
    ylims = (-75, -50),
    bottommargin = 10Plots.mm,
    topmargin = 3mm,
)
zz2 = plot(
    qq2,
    qq1,
    layout = (2, 1),
    size = (700, 400),
    legend = false,
    xlims=(-50,2050),
    leftmargin = 5Plots.mm,
)

layout = @layout [a{0.48w} _ b{0.48w}]
zz = plot(zz1, zz2, layout = layout, size = (1400, 400), frame=:axes)
path = mkpath(plotsdir("up-down", "Figures"))
savefig(zz, joinpath(path, "Fig7_experimental.pdf"))

zz

##
layout = @layout [
    a{0.7h}
    b{0.3h}
]
zzz = plot(pp, zz, layout = layout, size = (1400, 1100))

savefig(zzz, joinpath(path, "Fig7.pdf"))
zzz


# ##
# data = load(datadir("up_down","bistability","characterize_updown_(150, 150).jld2")) |> dict2ntuple
# @unpack βs, νs = data
# data.data
# yy1 = vcat([x.up[1,:] for x in data.data[1,5,5:3:end]],
# [x.up[1,:] for x in data.data[2,5,5:3:end]],
# [x.up[1,:] for x in data.data[3,5,5:3:end]])
# yy2 = vcat([x.down[1,:] for x in data.data[1,5,5:3:end]],
# [x.down[1,:] for x in data.data[2,5,5:3:end]],
# [x.down[1,:] for x in data.data[3,5,5:3:end]])
# xx = vcat(collect(3:3:44), collect(4:3:44 ), collect(5:3:44 ))
# c= repeat([:darkblue, :darkred, :black],inner=14)
# p = plot()
# plot!(p, xticks=(10:10:44, round.(νs[10:10:44],digits=0)), xlabel="Input rate (kHz)", ylabel="Av. duration (ms)", label="")
# plot!(p,yticks=(-300:150:300, ["300", "150", "0", "150", "300"]),label="")
# q = hline!(twinx(), [-50,50], ylims=(-300,300), lc=:black, ls=:dot, legend=false, yticks=([-50,50],[string(L"τ_r"),string(L"τ_r")] ),alpha=0.2, label="")
# plot!(twinx(), ylims=(-300,300), yticks=([-190,190], ["DOWN", "UP"]), legend=false, label="", yrotation=-90)
# bar!(p, xx,[mean.(yy1)], c=c, ylims=(0,10000), legend=false, lc=:white, label="")
# bar!(xx,-[mean.(yy2)], c=c, ylims=(-2000,2000), legend=false, lc=:white, label="")
# bar!([[],[],[]], c=c, legend=:bottomright, lc=:white, fg_color_legend=:transparent)
# p1 = plot!(size=(700,400), title="\nProximal - proximal")
# ##

# # analysis: length, rate, mem_avg, mem_std
# p1 = dotplot(1:6, [NMDA.up[1,:], NMDA.down[1,:], AMPA.up[1,:], AMPA.down[1,:], SOMA.up[1,:], SOMA.down[1,:] ], yscale=:log, ylabel="Length (ms)", xticks=([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN", "UP", "DOWN"]), legend=false, c = [:darkblue :darkblue :darkred :darkred :black :black], xrotation=-45, linewidth=0)
# # p1 = dotplot!(1:6, [NMDA.up[1,:], NMDA.down[1,:], AMPA.up[1,:], AMPA.down[1,:], SOMA.up[1,:], SOMA.down[1,:] ], yscale=:log, ylabel="Length (ms)", xticks=([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN", "UP", "DOWN"]), legend=false, c = [:darkblue :darkblue :darkred :darkred :black :black], xrotation=-45, linewidth=0)
# # p1 = scatter!([fill(length(NMDA.up[1,:]),1), length(NMDA.down[1,:]),  ],[NMDA.up[1,:], NMDA.down[1,:] ])
# # p1 = scatter!([1:length(NMDA.up[1,:]), 1:length(NMDA.down[1,:]),  ],[NMDA.up[1,:], NMDA.down[1,:], AMPA.up[1,:], AMPA.down[1,:], SOMA.up[1,:], SOMA.down[1,:]])
# # p1 = boxplot!(1:6,[NMDA.up[1,:], NMDA.down[1,:], AMPA.up[1,:], AMPA.down[1,:], SOMA.up[1,:], SOMA.down[1,:] ], yscale=:log, ylabel="Length (ms)", xticks=([1, 2, 3, 4, 5, 6], ["UP", "DOWN", "UP", "DOWN", "UP", "DOWN"]), legend=false, c = [:darkblue :darkblue :darkred :darkred :black :black], xrotation=-45, fillalpha=0.75)
# scatter!(mean.([NMDA.up[1,:], NMDA.down[1,:], AMPA.up[1,:], AMPA.down[1,:], SOMA.up[1,:], SOMA.down[1,:]]), yscale=:log, ylabel="Length (ms)", legend=false, c = :white, ms=6)
# # hline!([100], ylims=(10,1000), ls=:dash, c=:black, ms=4)
# ##  

# # dunction cumulative_probability()
# ##

# # states = NMDA.up[1,:] 
# # states = AMPA.up[1,:] 
# # states = SOMA.up[1,:] 
# # _hist = length.(UncertainData.bin(0:1:100_00,states  ,1:length(states)))
# # plot!(cum_prob(_hist), size=(300,400))
# # scatter(AMPA.up[2,nonzero],AMPA.up[4,nonzero],ylabel="Firing rate (Hz)", c=:darkblue, title="NMDA", msc=:darkblue, ms=4)
# # cor(NMDA.up[2,nonzero],NMDA.up[4,nonzero])


# data = load(datadir("up_down","bistability","characterize_updown_(400, 400).jld2")) |> dict2ntuple
# @unpack βs, νs = data
# data.data
# yy1 = vcat([x.up[1,:] for x in data.data[1,5,5:3:end]],
# [x.up[1,:] for x in data.data[2,5,5:3:end]],
# [x.up[1,:] for x in data.data[3,5,5:3:end]])
# yy2 = vcat([x.down[1,:] for x in data.data[1,5,5:3:end]],
# [x.down[1,:] for x in data.data[2,5,5:3:end]],
# [x.down[1,:] for x in data.data[3,5,5:3:end]])
# xx = vcat(collect(3:3:44), collect(4:3:44 ), collect(5:3:44 ))
# c= repeat([:darkblue, :darkred, :black],inner=14)
# p = plot()
# plot!(p, xticks=(10:10:44, round.(νs[10:10:44],digits=0)), xlabel="Input rate (kHz)", ylabel="Av. duration (ms)", label="")
# plot!(p,yticks=(-300:150:300, ["300", "150", "0", "150", "300"]),label="")
# q = hline!(twinx(), [-50,50], ylims=(-300,300), lc=:black, ls=:dot, legend=false, yticks=([-50,50],[string(L"τ_r"),string(L"τ_r")] ),alpha=0.2, label="")
# plot!(twinx(), ylims=(-300,300), yticks=([-190,190], ["DOWN", "UP"]), legend=false, label="", yrotation=-90)
# bar!(p, xx,[mean.(yy1)], c=c, ylims=(0,300), legend=false, lc=:white, label="")
# bar!(xx,-[mean.(yy2)], c=c, ylims=(-300,300), legend=false, lc=:white, label="")
# bar!([[],[],[]], c=c, legend=:bottomright, lc=:white, fg_color_legend=:transparent)
# p2 = plot!(size=(700,400), title="\nDistal - distal")
# pp3 = plot(p1,p2, layout=(1,2), size=(1400,400), legend=false, leftmargin=5Plots.mm,)

# layout= @layout [a{0.5w}	 b{0.5w}
# 				a{0.3h}]
# zz =plot(pp1,pp2,pp3, layout=layout, size=(1400,1000), )
# ##

# path = mkpath(plotsdir("up_down","Figures"))
# zz
# savefig(zz, joinpath(path,"Fig7.pdf"))
