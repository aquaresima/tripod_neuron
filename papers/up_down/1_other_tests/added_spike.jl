##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
include(projectdir("scripts", "plots", "default_plots.jl"))


##
get_data() = load(datadir("up_down", "activity", "added_spike.jld2"))
added_spike = (; Dict(Symbol(k) => x for (k, x) in get_data())...)
data = @unpack spiketime, dist_NMDA, prox_NMDA, dist_AMPA, prox_AMPA, soma_AMPA, soma_NMDA =
    added_spike

@unpack νs, ds = TN


##
# Plots Spike triggered averages
default(palette = palette(:tab10))
c = collect(palette(:roma, size(νs[TN.min_NMDA:end])...))'
layout = @layout[a{0.44w} b{0.44w} _]
q11 = plot(
    prox_NMDA.mem[1, spiketime-50:end, TN.min_NMDA:end],
    c = c,
    legend = false,
    ylabel = " ΔEPSP (mV)",
    title = "proximal",
    xlabel = " Time (ms)",
)

# scatter!(prox_NMDA.half[TN.min_NMDA:end] ..-500, ones(28), c=c')
q12 = plot(
    dist_NMDA.mem[1, spiketime-50:end, TN.min_NMDA:end],
    c = c,
    legend = false,
    xlabel = "Time (ms)",
    ylabel = " ΔEPSP (mV)",
    title = "distal",
)


q1 = plot(
    q11,
    q12,
    xticks = (50:1050:3050, 0:100:300),
    ylims = (0, 1),
    xlims = (0, 3000),
    frame = :origin,
    layout = layout,
)
hm = zeros(length(TN.νs[TN.min_NMDA:end]), 1) |> x -> (x[:, 1] .= 1:length(x); x)
inset = (1, bbox(0.0, 0.4, 0.1, 0.30, :bottom, :right))
xs = round.(Int, range(1, 28, 3))
xxs = round.(Int, TN.νs[collect(xs .+ TN.min_NMDA .- 1)])
_xticks = (xs, xxs)
q1 = heatmap!(
    hm,
    c = :roma,
    frame = :box,
    xticks = :none,
    cbar = false,
    inset = inset,
    subplot = 3,
    yticks = _xticks,
    ylabel = "Input rate (kHz)",
    guidefontsize = 15,
    tickfontsize = 15,
)


ms = 10
q21 = scatter(
    TN.νs[TN.min_NMDA:end-7],
    soma_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = :black,
    msc = :black,
    shape = :square,
    ms = 10,
    label = "soma",
)
q21 = scatter!(
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = nmda_color,
    msc = nmda_color,
    ms = ms,
    label = "proximal NMDA",
)
q21 = scatter!(
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = ampa_color,
    msc = ampa_color,
    ms = ms,
    label = "proximal AMPA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = nmda_color,
    msc = nmda_color,
    shape = :diamond,
    ms = ms,
    label = "distal NMDA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = ampa_color,
    msc = ampa_color,
    shape = :diamond,
    ms = ms,
    label = "distal AMPA",
)
q2 = plot!(ylabel = "Width (ms)", xlabel = "Input rate (kHz)")
plot!(legend = false, ylims = (0, 80))
qq = plot()

scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    soma_AMPA.peak[TN.min_NMDA:end-7],
    c = :transparent,
    msc = :black,
    ms = 10,
    shape = :square,
    label = "Soma only",
    yticks = (0:0.5:1.5, 0:0.5:1.5),
    ylabel = "Peak (mV)",
    ylims = (0, 1.5),
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.peak[TN.min_NMDA:end-7],
    c = nmda_color,
    msc = nmda_color,
    ms = ms,
    label = "proximal NMDA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.peak[TN.min_NMDA:end-7],
    c = nmda_color,
    msc = nmda_color,
    shape = :diamond,
    ms = ms,
    label = "distal NMDA",
    yticks = (0:1.0, 0:1),
    ylabel = "Peak (mV)",
    ylims = (0, 1),
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.peak[TN.min_NMDA:end-7],
    c = ampa_color,
    msc = ampa_color,
    ms = ms,
    label = "proximal AMPA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.peak[TN.min_NMDA:end-7],
    c = ampa_color,
    msc = ampa_color,
    ms = ms,
    shape = :diamond,
    label = "distal AMPA",
    ylabel = "Peak (mV)",
    ylims = (0, 1.5),
    xlabel = "Input rate (kHz)",
)
plot!(legendfontsize = 13, legend = :topright)


# leg_colors = [ nmda_color nmda_color ampa_color   ampa_color :black]
# # q4 = scatter(frame=:none,[[1],[1],[1],[1],[1]],labels=["  prox NMDA" "  dist NMDA" "  prox AMPA"  "  dist AMPA" "  soma only"], c= leg_colors, msc= leg_colors, ms=[3 3 3 3 2], shape=[:circle  :diamond :circle :diamond :square], ylims=(2,2))
# # plot!(legendfontsize=15,)




layout = @layout[a{0.45w} b{0.45w}]
q_2 = plot(q2, qq, layout = layout, xscale = :log, ticksfontsize = 15)
q = plot(
    q1,
    q_2,
    layout = (2, 1),
    size = two_columns_size,
    leftmargin = 10Plots.mm,
    bottommargin = 10Plots.mm,
)
##
