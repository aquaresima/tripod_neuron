##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars

include(projectdir("scripts", "plots", "default_plots.jl"))

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
ampa_color = :darkred
nmda_color = :darkblue

##
using Random
simtime = 3000
dd = (150, 150)
istdp_syn = false
do_spikes = true
_rate = 20
nmda = false
νs

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

layout = @layout [a b; c d; e f]
p = plot(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    lw = 2,
    legendfontsize = 13,
    layout = layout,
    legend_background_color = :transparent,
    size = two_columns_size .* (1.0, 1.5),
)

##
new_ampa = my_diag(AMPA.rate)
new_ampa[3, 5] ≈ AMPA.rate[1, 3, 3, 5]
my_diag(AMPA.std_v, 2)


## Average activity
default(margin = 0Plots.mm)
soma_rate = ((AMPA.rate[:, 1:1, 1:1, :].+NMDA.rate[:, 1:1, 1:1, :])./2)[:, 1, 1, :]
q_bis = plot(
    TN.νs,
    soma_rate[1, :],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "soma",
    c = :black,
    ylims = (0, 2),
)
q_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate[:, 2:end, 2:end, :])', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate)', dims = 2),
    label = "NMDA",
    c = nmda_color,
)
q_bis = plot!(
    TN.νs,
    mean(my_diag(AMPA.rate[:, 2:end, 2:end, :])', dims = 2)[:, 1],
    ribbon = std(my_diag(AMPA.rate)', dims = 2),
    label = "AMPA",
    c = ampa_color,
)
plot!(xscale = :log, yticks = 0:2:6)
plot!(ylims = (-2, 6), ylabel = "Firing \nrate (Hz)", legend = false, xticks = :false)
plot(yticks = :none)
plot!(ylims = (-1, 2.5), ylabel = "ISI CV", legend = false, xticks = :false)
q2_bis = plot!(
    TN.νs,
    soma_rate[2, :],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "soma",
    c = :black,
)
q2_bis = plot!(
    TN.νs,
    mean(my_diag(AMPA.rate[:, 2:end, 2:end, :], 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "AMPA",
    c = ampa_color,
)
q2_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate[:, 2:end, 2:end, :], 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "NMDA",
    c = nmda_color,
)
plot!(xscale = :log, legend = false, yticks = 0:1:2)
# q78 = plot(p_bis, p2_bis, xscale=:log,link=:x, lw=3, layout=(1,2), xlabel="Input rates (kHz)", leftmargin=10Plots.mm)

q3_bis = plot(
    νs,
    mean(NMDA.std_v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    ribbon = std(NMDA.std_v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    c = :darkblue,
)
plot!(
    νs,
    mean(AMPA.std_v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    ribbon = std(AMPA.std_v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :],
    c = :darkred,
)
plot!(νs, AMPA.std_v[1, 1, 1, :], c = :black, legend = false)
plot!(xscale = :log)
plot!(xlabel = "Input rates (kHz)")
vline!([νs[20]], c = :black, linestyle = :dash, label = "")
q_l = plot(
    [[], [], []],
    frame = :none,
    c = [:black ampa_color nmda_color],
    label = ["  Soma" "  AMPA only" "  NMDA"],
    legendfontsize = 18,
    legend = :left,
    legend_background_color = :transparent,
)
layout = @layout [grid(3, 1) c{0.2w} _]
q = plot(q_bis, q2_bis, q3_bis, q_l, layout = layout, size = two_columns_size .* (0.7, 1.5))
plot!(bottommargin = 5Plots.mm, topmargin = 5Plots.mm)

# savefig(zD, joinpath(path,"Fig4D.pdf"))
# zD
##
##


s1 = heatmap(
    TN.νs,
    TN.ds[2:end],
    my_diag(NMDA.rate)[2:end, :],
    title = "AMPA + NMDA",
    clims = (0, 10),
    c = :blues,
    cbar = false,
    xscale = :log,
    xticks = :none,
    yticks = [100, 300, 500],
)
# s2 = heatmap(TN.νs, TN.ds[2:end], my_diag(AMPA.rate)[2:end,:],title="AMPA only", clims=(0,10), c=:blues, cbar=false, xscale=:log, yticks=:none, xticks=:none)
s4 = heatmap(
    TN.νs,
    TN.ds[2:end],
    my_diag(NMDA.std_v, 1)[2:end, :],
    clims = (0, 6),
    cbar = false,
    xscale = :log,
    xlabel = "                   Input rate (kHz)",
    yticks = [100, 300, 500],
)
# s5 = heatmap(TN.νs, TN.ds[2:end], my_diag(AMPA.std_v)[2:end,:],clims=(0,6),yticks=:none, xscale=:log, cbar=false)
# _s4 = twinx(s5)
plot!(s4, ylabel = "                             Dendritic length (μm)")


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
layout = @layout[[a{0.480w} aa{0.05w}]; [a{0.46w} aa{0.05w}]]
s = plot(
    s1,
    ss1,
    s4,
    ss2,
    colorbar_tickfontsize = 11,
    layout = layout,
    size = (700, 600 * 7 / 9),
)

# plot!(bottommargin=5Plots.mm, topmargin=5Plots.mm)
##

layout = @layout [
    a{0.38h}
    _
    b{0.58h}
    _
]
qs = plot(q, s, layout = layout, size = two_columns_size .* (1, 1))
l = @layout [a{0.45w} _ b{0.40w,1.04h} _]
# ]
z = plot(p, qs, layout = l, size = two_columns_size .* (1.2, 1.2), margin = 5Plots.mm)

z
path = mkpath(plotsdir("up_down", "Figures"))
savefig(z, joinpath(path, "Fig4.pdf"))
z


##


simtime = 1000
function run_conductance_ex(nmda, dd; ld = true, kwargs...)
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
    voltage, syn = TN.run_tripod(
        dd,
        inputs,
        simtime;
        synapses = cond.synapses,
        model...,
        rec_syn = true,
    )
    ld1 = ld ? "Dendrite $(dd[1]) μm" : ""
    ld2 = ld ? "Dendrite $(dd[2]) μm" : ""
    ls = ld ? "" : "Soma"
    voltage[voltage.>-40] .= -40
    scale = 1e-3
    syn = -syn
    p = plot(
        scale .* map(
            t -> TN.syn_current_glu(voltage[2, t], syn[2, t, :], syn_model.Esyn_dend),
            1:simtime*10,
        ),
        palette = :default,
        label = "Glu",
    )
    plot!(
        scale .* map(
            t -> TN.syn_current_gaba(voltage[2, t], syn[2, t, :], syn_model.Esyn_dend),
            1:simtime*10,
        ),
        palette = :default,
        label = "Gaba",
    )
    plot!(
        scale .* map(
            t -> TN.syn_current_tot(voltage[2, t], syn[2, t, :], syn_model.Esyn_dend),
            1:simtime*10,
        ),
        palette = :default,
        c = :black,
        ls = :dash,
        label = "Total",
    )
    plot!(ylims = (-4, 4), xlims = (simtime / 3 * 2 * 10, simtime * 10))
    plot!(xticks = (20_000:2500:30000, 0:0.25:1), yticks = ([-3, 0, 3]))
    plot!(; kwargs...)
    return p
end

##
p1 = run_conductance_ex(true, (200, 200), legend = false)
p2 = run_conductance_ex(false, (200, 200), yticks = :none)
z = plot(p1, p2, size = two_columns_size .* (0.5, 0.3), margin = 5Plots.mm)
savefig(z, joinpath(path, "Fig4E.pdf"))
