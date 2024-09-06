##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
include(projectdir("scripts", "up_down", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_plots.jl"))
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
    plot!(ylims = (-80, 80))
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
p5 = run_example(
    true,
    (400, 400),
    xlabel = "                                 Time (s)",
    ld = false,
    xticks = (0:2500:10000, 0:0.25:1),
    yticks = (-70:25:25, -70:25:25),
)
p6 = run_example(false, (400, 400), yticks = :none, xticks = (0:2500:10000, 0:0.25:1))

p = plot(
    p1,
    p2,
    p5,
    p6,
    lw = 2,
    size = (600, 500),
    layout = layout,
    legend_background_color = :transparent,
)

path = mkpath(plotsdir("up_down", "poster"))
savefig(p, joinpath(path, "dendritic_dynamics.pdf"))


##
p_bis = plot(
    TN.νs,
    mean(my_diag(AMPA.rate)', dims = 2)[:, 1],
    ribbon = std(my_diag(AMPA.rate)', dims = 2),
    label = "AMPA",
    c = :darkred,
    xlabel = "                                  Input rate (kHz)",
)
p_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate)', dims = 2),
    label = "NMDA",
    ylims = (0, 20),
    ylabel = "Firing rate (Hz)",
    legend = :topleft,
    c = :darkblue,
)
p2_bis = plot(
    TN.νs,
    mean(my_diag(AMPA.rate, 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "AMPA",
    c = :darkred,
    ylims = (0, 2),
)
p2_bis = plot!(
    TN.νs,
    mean(my_diag(NMDA.rate, 2)', dims = 2)[:, 1],
    ribbon = std(my_diag(NMDA.rate, 2)'),
    label = "NMDA",
    c = :darkblue,
    ylims = (0, 1.5),
    ylabel = "ISI CV",
    legend = false,
)
layout = @layout [grid(2, 2)]
default(size = (600, 300), margin = 5mm)
p78 = plot(p_bis, p2_bis, xscale = :log, lw = 3, legend = false, xlabel = "")
##

q21 = scatter(
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = :darkblue,
    msc = :darkblue,
    ms = 5,
    label = "proximal NMDA",
)
q21 = scatter!(
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = :darkred,
    msc = :darkred,
    ms = 5,
    label = "proximal AMPA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.half[TN.min_NMDA:end-7] .- 500,
    c = :darkblue,
    msc = :darkblue,
    shape = :star,
    ms = 8,
    label = "distal NMDA",
)
scatter!(
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.half[TN.min_NMDA:end-7] .- 500,
    c = :darkred,
    msc = :darkred,
    shape = :star,
    ms = 8,
    label = "distal AMPA",
)
q2 = plot!(ylabel = "Width (ms)", xlabel = "                             Input rate (kHz)")
plot!(legend = false)
#
#
qq = plot()
# qq = twinx(q3)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_NMDA.peak[TN.min_NMDA:end-7],
    c = :darkblue,
    msc = :darkblue,
    ms = 5,
    label = "proximal NMDA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    prox_AMPA.peak[TN.min_NMDA:end-7],
    c = :darkred,
    msc = :darkred,
    ms = 5,
    label = "proximal AMPA",
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_NMDA.peak[TN.min_NMDA:end-7],
    c = :darkblue,
    msc = :darkblue,
    shape = :star,
    ms = 5,
    label = "distal NMDA",
    yticks = (0:1.0, 0:1),
    ylabel = "Peak (mV)",
    ylims = (0, 1),
)
scatter!(
    qq,
    TN.νs[TN.min_NMDA:end-7],
    dist_AMPA.peak[TN.min_NMDA:end-7],
    c = :darkred,
    msc = :darkred,
    ms = 5,
    shape = :star,
    label = "distal AMPA",
    ylabel = "Peak (mV)",
    ylims = (0, 1),
)
# scatter!(qq, TN.νs[TN.min_NMDA:end-7],soma_AMPA.peak[TN.min_NMDA:end-7], c=:black, msc=:black, ms=3, shape=:square,label="distal AMPA", yticks=(0:1., 0:1), ylabel="Peak (mV)", ylims=(0,1))
plot!(legendfontsize = 4)

# q4 = scatter(frame=:none,[[1],[1],[1],[1],[1]],labels=["  prox NMDA" "  dist NMDA" "prox AMPA"  "dist AMPA" "soma only"], c=[:darkred   :darkred  :darkblue :darkblue :black], msc=[:darkred :darkred :darkblue :darkblue :black], ms=[3 3 3 3 2], shape=[:circle :circle :star :star :square], ylims=(2,2))
# plot!(legendfontsize=10)




q_2 = plot(q2, qq, xscale = :log, link = :x, ticksfontsize = 9, size = (600, 300))
p = plot(p78, q_2, layout = (2, 1), size = (600, 500))

path = mkpath(plotsdir("up_down", "poster"))
savefig(p, joinpath(path, "firing_EPSP.pdf"))
##

##

function get_trace(model, β, rate; do_spikes = false)
    simtime = 25_000
    seed = 39
    AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)

    cond = TN.get_balance_conditions(model.ds..., rate; nmda = model.nmda)
    inputs = TN.make_spikes(simtime, β; cond..., seed = seed)
    voltage = TN.run_tripod(
        inputs,
        simtime;
        model...,
        do_spikes = true,
        syn_model = cond.syn_model,
        AdEx = AdEx,
    )
    pp = plot(voltage[1, end-20000:end], c = :gray)
    if model.nmda
        c = :darkblue
        if β == 50
            plot!(pp, title = "AMPA + NMDA")
        end
    else
        c = :darkred
        if β == 50
            plot!(pp, title = "AMPA only")
        end
    end
    plot!(voltage[2, end-20000:end], c = c)
    plot!(
        ylims = (-60, 0),
        xlabel = "Time (s)",
        yticks = (-60:30:0, -60:30:0),
        frame = :axes,
        xticks = (0:10_000:20_000, 0:1:2),
        ylabel = "                        Membrane potential (mV)",
        legend = false,
    )
    aph = APH(voltage, "", filter = false)
    plot!(
        aph,
        frame = :axes,
        xlabel = "Soma (mV)",
        xrotation = -45,
        title = "β = $β",
        ylabel = "                             Bin frequency ",
    )
    plot(pp, aph, layout = (1, 2), right_margin = 5Plots.mm)


    # TN.get_balance_conditions(model.ds...,rate, nmda=model.nmda)
    # inputs = TN.make_spikes(simtime, β; cond...)
    # voltage = TN.run_tripod(inputs,simtime, do_spikes=do_spikes; model..., syn_model=cond.syn_model)
    # output = round(TN.get_spike_rate(voltage[1,:]), digits=2)
    # μ = mean(voltage[:,5000:end], dims=2)[:,1]
    # plot(voltage[1,:])

    # p1 = plot!(xticks=:none)
    # p1 = plot!(subplot=3,xticks=(range(1, simtime*10, length=5), round.(range(1, simtime/1000, length=5),digits=2)), xlabel="Time (s)")
end
# plot!(ylims=(50_000, 70_000))

#
rate = 18
model = (ds = (300, 300), nmda = true)
NMDAplots = []
for β in [50, 150]
    p = get_trace(model, β, rate, do_spikes = true)
    push!(NMDAplots, p)
end
model = (ds = (300, 300), nmda = false)
AMPAplots = []
for β in [50, 150]
    p = get_trace(model, β, rate, do_spikes = true)
    push!(AMPAplots, p)
end

# TN.ampa_equivalent.Esyn_dend.AMPA.g0
for (p, q) in zip(NMDAplots[1:end-1], AMPAplots[1:end-1])
    plot!(p, frame = :axes, xlabel = "", xtick = :none, ylabel = " ")
    plot!(q, frame = :axes, xlabel = "", xtick = :none, ylabel = " ")
end

layout = @layout[
    a{0.5h}
    a{0.5h}
]
p = plot(
    plot(NMDAplots[1:2]..., layout = layout, legend = false, size = (1000, 500)),
    plot(AMPAplots[1:2]..., layout = layout, legend = false, size = (1000, 500)),
    legend = false,
    tickfontsize = 11,
)

path = mkpath(plotsdir("up_down", "poster"))
savefig(p, joinpath(path, "fluctuations.pdf"))


##
##
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
simtime = 25_000
rate = 18;
β = 240;
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
    titlefontsize = 15,
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
        label = "$w mV",
        legendtitle = "Window",
        legend_fontsize = 19,
        legend = :topleft,
    )
end

m = maximum(globalKDE(1, voltage[1, :], v_range = -80:-40))
plot!(yticks = :none, frame = :axes, yaxis = false, margin = 2Plots.mm, ylims = (0, m))

hm = zeros(50, 1)
hm[:, 1] .= collect(1:length(hm))
inset = (1, bbox(0.0, 0.60, 0.1, 0.3, :bottom, :right))
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

get_data(syn) = load(datadir("up_down", "robustness", "critical_window_$syn.jld2"))

AMPA = get_data("AMPA")["data"]
NMDA = get_data("NMDA")["data"]
βs = get_data("NMDA")["βs"]

titles = ["distal\ndistal", "distal\nproximal", "proximal\nproximal", "soma\nonly"]

plotsNMDA = []
for (n, (k, v)) in enumerate(NMDA)
    p = heatmap(TN.νs, βs, v.critical', clims = (0, 50), c = :amp)
    push!(plotsNMDA, p)
    m = argmax(mean(v.critical[:, :], dims = 2))[1]
    plot!(title = titles[n])
    # vline!([TN.νs[m]])
    @show m
end

plotsAMPA = []
for (n, (k, v)) in enumerate(AMPA)
    p = heatmap(TN.νs, βs, v.critical', clims = (0, 50), c = :amp)
    push!(plotsAMPA, p)
    m = argmax(mean(v.critical[:, :], dims = 2))[1]
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
##
# layout = @layout [a{0.95h, 0.35w} _ a{0.6w} _]
p = plot(z1, z2, size = (600, 1000), layout = (2, 1), topmargin = 10Plots.mm)

path = mkpath(plotsdir("up_down", "poster"))
savefig(p, joinpath(path, "up_down.pdf"))
##
