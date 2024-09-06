##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
using EllipsisNotation

include(projectdir("scripts", "plots", "default_plots.jl"));

##
get_data(title, type = "data") =
    load(datadir("up_down", "activity", "full_activity_$title.jld2"))
AMPA = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)

get_data() = load(datadir("up_down", "activity", "added_spike.jld2"))
added_spike = (; Dict(Symbol(k) => x for (k, x) in get_data())...)
data = @unpack spiketime, dist_NMDA, prox_NMDA, dist_AMPA, prox_AMPA, soma_AMPA, soma_NMDA =
    added_spike

@unpack νs, ds = TN

##
using Random
simtime = 3000
dd = (150, 150)
istdp_syn = false
do_spikes = true
_rate = 20
nmda = false

function run_example(nmda, dd; ld = true, kwargs...)
    Random.seed!(17)
    syn_model = nmda ? TN.human_synapses : TN.ampa_equivalent
    model = (
        ds = dd,
        syn_model = syn_model,
        species = "H",
        soma_only = dd[1] == 0 ? true : false,
        do_spikes = do_spikes,
    )
    cond = TN.get_balance_conditions(dd..., _rate, istdp_syn = istdp_syn, nmda = nmda)
    inputs = TN.make_spikes(simtime; cond...)
    sum(inputs, dims = 2)
    plot(; kwargs...)
    voltage = TN.run_tripod(dd, inputs, simtime; synapses = cond.synapses, model...)
    ld1 = ld ? "Dendrite $(dd[1]) μm" : ""
    ld2 = ld ? "Dendrite $(dd[2]) μm" : ""
    ls = ld ? "" : "Soma"
    _p = plot!(voltage[2, end-10_000+1:end], label = ld1)
    _p = plot!(voltage[3, end-10_000+1:end], label = ld2)
    _p = plot!(voltage[1, end-10_000+1:end], c = :black, label = ls)
    plot!(ylims = (-80, 80), legend = :topright)
    # p=	plot!(xlims=(2e4,3e4))
    return _p
end
p1 = run_example(
    true,
    (150, 150),
    xticks = :none,
    title = "AMPA + NMDA",
    ld = false,
    yticks = (-70:25:25, -70:25:25),
)
p2 = run_example(false, (150, 150), xticks = :none, title = "AMPA only", yticks = :none)
p3 = run_example(
    true,
    (400, 400),
    xticks = :none,
    ylabel = "Membrane potential (mV)",
    ld = false,
    yticks = (-70:25:25, -70:25:25),
)
p4 = run_example(false, (400, 400), xticks = :none, yticks = :none)
p5 = run_example(
    true,
    (150, 400),
    xlabel = "Time (s)",
    ld = false,
    xticks = (0:2500:10000, 0:0.25:1),
    yticks = (-70:25:25, -70:25:25),
)
p6 = run_example(
    false,
    (150, 400),
    yticks = :none,
    xlabel = "Time (s)",
    xticks = (0:2500:10000, 0:0.25:1),
)

AMPA.rate
soma_rate = ((AMPA.rate[:, 1:1, 1:1, :].+NMDA.rate[:, 1:1, 1:1, :])./2)[:, 1, 1, :]
p_bis = plot(
    TN.νs,
    soma_rate[1, :],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "soma",
    c = :black,
    ylims = (0, 2),
)
p_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate[:, 2:end, 2:end, :])', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate)', dims = 2),
    label = "NMDA",
    ylims = (0, 6),
    ylabel = "Firing rate (Hz)",
    legend = :topleft,
    c = nmda_color,
)
p_bis = plot!(
    TN.νs,
    mean(my_diag(AMPA.rate[:, 2:end, 2:end, :])', dims = 2)[:, 1],
    ribbon = std(my_diag(AMPA.rate)', dims = 2),
    label = "AMPA",
    c = ampa_color,
    xlabel = "                                  Input rates (kHz)",
)
p2_bis = plot(
    TN.νs,
    soma_rate[2, :],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "soma",
    c = :black,
    ylims = (0, 2),
)
p2_bis = plot!(
    TN.νs,
    mean(my_diag(AMPA.rate[:, 2:end, 2:end, :], 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "AMPA",
    c = ampa_color,
    ylims = (0, 2),
)
p2_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate[:, 2:end, 2:end, :], 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "NMDA",
    c = nmda_color,
    ylims = (0, 2),
    ylabel = "ISI CV",
    legend = false,
)
p78 = plot(
    p_bis,
    p2_bis,
    xscale = :log,
    link = :x,
    lw = 3,
    layout = (1, 2),
    xlabel = "Input rates (kHz)",
    leftmargin = 10Plots.mm,
)

layout = @layout [
    grid(3, 2)
    a{0.15h}
]

p = plot(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p78,
    lw = 2,
    legendfontsize = 13,
    layout = layout,
    legend_background_color = :transparent,
    size = two_columns_size .* (1.0, 1.5),
)
##
# Plots Spike triggered averages
default(palette = palette(:tab10), margins = 2mm)
c = collect(palette(:roma, size(νs[TN.min_NMDA:end])...))'
q11 = plot(
    prox_NMDA.mem[1, spiketime-50:end, TN.min_NMDA:end],
    c = c,
    legend = false,
    ylabel = " ΔEPSP (mV)",
    title = "proximal",
    xlabel = "                             Time (ms)",
)
plot!(xticks = (50:1050:3050, 0:100:300), ylims = (0, 0.8), xlims = (0, 3000))
plot(ticks = :none)
q12 = plot!(
    twinx(),
    dist_NMDA.mem[1, spiketime-50:end, TN.min_NMDA:end],
    c = c,
    legend = false,
    title = "distal",
)
plot!(xticks = (50:1050:3050, 0:100:300), ylims = (0, 0.8), xlims = (0, 3000))
q1 = colorbar(νs[1], νs[end], 4, c = :roma, title = "Input rate (kHz)")
layout = @layout[a{0.45w} c{0.1w} _ b{0.45w}]
q1 = plot(q11, q1, q12, layout = layout)
##
q21 = scatter(
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = nmda_color,
    msc = nmda_color,
    ms = 5,
    label = "proximal NMDA",
)
q21 = scatter!(
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = ampa_color,
    msc = ampa_color,
    ms = 5,
    label = "proximal AMPA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = nmda_color,
    msc = nmda_color,
    shape = :star,
    ms = 8,
    label = "distal NMDA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = ampa_color,
    msc = ampa_color,
    shape = :star,
    ms = 8,
    label = "distal AMPA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    soma_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = :black,
    msc = :black,
    shape = :square,
    ms = 3,
    label = "soma",
)
q2 = plot!(
    ylabel = "Width (ms)",
    yticks = (10:20:70, 10:20:70),
    xlabel = "                             Input rate (kHz)",
)
plot!(legend = false)

plot(ticks = :none)
qq = plot!(
    twinx(),
    yticks = (0:0.5:1.5, 0:0.5:1.5),
    ylims = (0, 1.5),
    frame = :box,
    ylabel = "Peak (mV)",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.peak[TN.min_NMDA:end-7],
    c = nmda_color,
    msc = nmda_color,
    ms = 5,
    label = "proximal NMDA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.peak[TN.min_NMDA:end-7],
    c = ampa_color,
    msc = ampa_color,
    ms = 5,
    label = "proximal AMPA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.peak[TN.min_NMDA:end-7],
    c = nmda_color,
    msc = nmda_color,
    shape = :star,
    ms = 5,
    label = "distal NMDA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.peak[TN.min_NMDA:end-7],
    c = ampa_color,
    msc = ampa_color,
    ms = 5,
    shape = :star,
    label = "distal AMPA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    soma_AMPA.peak[TN.min_NMDA:end-7],
    c = :black,
    msc = :black,
    ms = 3,
    shape = :square,
    label = "distal AMPA",
)
plot!(legendfontsize = 4, legend = false)

q4 = scatter(
    frame = :none,
    [[1], [1], [1], [1], [1]],
    labels = ["  prox NMDA" "  dist NMDA" " prox AMPA" " dist AMPA" " soma only"],
    c = [ampa_color ampa_color nmda_color nmda_color :black],
    msc = [ampa_color ampa_color nmda_color nmda_color :black],
    ms = [3 3 3 3 2],
    shape = [:circle :circle :star :star :square],
    ylims = (2, 2),
)
plot!(legendfontsize = 14)




layout = @layout[a{0.45w} c{0.2w} b{0.45w}]
q_2 = plot(q2, q4, qq, layout = layout, xscale = :log, link = :x, ticksfontsize = 9)
q = plot(q1, q_2, layout = (2, 1), size = (800, 600))
##


new_ampa = my_diag(AMPA.rate)
new_ampa[3, 5] ≈ AMPA.rate[1, 3, 3, 5]
my_diag(AMPA.std_v, 2)




s1 = heatmap(
    TN.ds[2:end],
    TN.νs,
    my_diag(NMDA.rate)[2:end, :]',
    title = "AMPA + NMDA",
    clims = (0, 10),
    c = :blues,
    cbar = false,
    yscale = :log,
    xticks = :none,
)
s2 = heatmap(
    TN.ds[2:end],
    TN.νs,
    my_diag(AMPA.rate)[2:end, :]',
    title = "AMPA only",
    clims = (0, 10),
    c = :blues,
    cbartitle = "output rate (Hz)",
    cbar = false,
    yscale = :log,
    yticks = :none,
    xticks = :none,
)
s4 = heatmap(
    TN.ds[2:end],
    TN.νs,
    my_diag(NMDA.std_v, 1)[2:end, :]',
    clims = (0, 6),
    cbar = false,
    yscale = :log,
    xlabel = "                   Dendritic length (μm)",
    xticks = (100:200:500, 100:200:500),
)
s5 = heatmap(
    TN.ds[2:end],
    TN.νs,
    my_diag(AMPA.std_v)[2:end, :]',
    clims = (0, 6),
    yticks = :none,
    yscale = :log,
    cbar = false,
    xticks = (100:200:500, 100:200:500),
)
# _s4 = twinx(s5)
plot!(s4, ylabel = "                             Input rate (kHz)", xrotation = 45)
plot!(s5, xrotation = 45)


hm = zeros(60, 1)
hm[:, 1] .= collect(1:length(hm))
inset = (1, bbox(-0.9, 1.1, 1, 1, :bottom, :right))
xs = round.(Int, range(1, 23, 4))
xxs = round.(Int, TN.νs[collect(xs .+ TN.min_NMDA .- 1)])
_xticks = (xs, xxs)
plot(frame = :none)
ss1 = heatmap!(
    twinx(),
    hm[1:10, :],
    c = :blues,
    cbar = false,
    xticks = :none,
    frame = :box,
    ylabel = "Output rate (Hz)",
    guidefontsize = 15,
)
plot(frame = :none)
ss2 = heatmap!(
    twinx(),
    [1],
    1:60,
    hm[1:1:end, :],
    cbar = false,
    xticks = :none,
    frame = :box,
    ylabel = "Membrane std (mV)",
    guidefontsize = 15,
    yticks = (0:10:60, 0:1:6),
)
#
s2
layout = @layout[[a{0.460w} b{0.46w} aa{0.05w} _]; [a{0.46w} b{0.46w} aa{0.05w} _]]
s = plot(
    s1,
    s2,
    ss1,
    s4,
    s5,
    ss2,
    colorbar_tickfontsize = 11,
    layout = layout,
    size = (700, 600 * 7 / 9),
)
##

layout = @layout [
    a{0.48h}
    _
    b{0.48h}
]
qs = plot(q, s, layout = layout, size = two_columns_size .* (1, 1.5))
l = @layout [a{0.45w} _ b{0.45w} _]
# ]
z = plot(p, qs, layout = l, size = two_columns_size .* (1.2, 1.5), margin = 5Plots.mm)



# path = mkpath(plotsdir("up_down","activity"))
# savefig(z, joinpath(path,"dendritic_dynamics.pdf"))
z

##
scatter(AMPA.mean_v[1, 2:end, 2:end, TN.min_AMPA:end][:])
mean(AMPA.mean_v[1, 2:end, 2:end, TN.min_AMPA:end][:])
mean(NMDA.mean_v[1, 2:end, 2:end, TN.min_NMDA:end][:])
scatter(NMDA.mean_v[1, 2:end-1, 2:end-1, TN.min_NMDA:end][1:9:end])
scatter(AMPA.mean_v[1, 2:end-1, 2:end-1, TN.min_AMPA:end][1:9:end])

plot(mean(AMPA.rate[1, 1, 1:1, :], dims = 1)[1, :])


heatmap(AMPA.mean_v[1, :, :, 24], clims = (-55, -45), rightmargin = 10Plots.mm)

plot(mean(NMDA.mean_v[1, :, :, :], dims = [1, 2])[1, 1, :])
plot!(mean(AMPA.mean_v[1, :, :, :], dims = [1, 2])[1, 1, :])
