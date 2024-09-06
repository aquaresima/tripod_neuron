using DrWatson
@quickactivate "Tripod"
include(projectdir("scripts", "up_down", "bimodal_kernel.jl"))
include(projectdir("scripts", "up_down", "updown_plots.jl"))
using TripodNeuron



##
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
simtime = 20_000
rate = 17
β = 200;
models_hist = []
nmda = true
ratio = 1.0
for (model, title) in zip(TN.models, TN.labels)
    cond = TN.get_balance_conditions(model.ds..., rate; nmda = true, shift = ratio)
    inputs = TN.make_spikes(simtime, β; cond...)
    voltage =
        TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)

    c = critical_window(voltage[1, :])
    aph = APH(voltage, title * " c: $c", 4, yticks = :none)
    r = TN.get_spike_rate(voltage[1, :])
    @show r, mean(voltage[1, :])
    push!(models_hist, aph)
end
p = plot(
    models_hist...,
    yaxis = false,
    tickfontsize = 13,
    rotation = -45,
    guidefontsize = 18,
)


z = plot(p, xrotation = -45, yticks = :none)
# savefig(z,joinpath(@__DIR__,"figures","all_models.pdf") )

##
kernel = globalKDE(16, voltage[1, :], collect(-70:-35))
plot(kernel)
isbimodal(kernel, 0.3)

voltage[1, :]
#         if !isbimodal(kernel, ratio)
#             return h
#         end
#     end
#     return max_b
# end
##
models_hist = []
β = 0.1
simtime = 100000
_rate = 1000

for prox = 0:2:6
    model = (150, 350)
    eiratio = [0.0, 0.5, 1.0]
    inputs = TN.make_spikes(
        simtime,
        β,
        baserate = 10.0,
        rate = _rate,
        soma = false,
        dends = true,
        gsyn = 1,
        eiratio = eiratio,
    )

    voltage =
        run_tripod(model, inputs, simtime, synapses = [0.0, prox, 1.0, 0.0, prox, 1.0])

    p = histogram(
        voltage[1, :],
        normalize = true,
        bins = range(-90, stop = -30, length = 180),
        color = :black,
        label = "",
        title = "",
        titlefontsize = 15,
    )
    s = maximum(map(x -> !isnan(x) ? x : 0, p.series_list[1][:y]))
    kde = globalKDE(3, voltage[1, :])
    plot!(-90:-40, kde, label = "", lw = 2)
    rate = round(TN.get_spike_rate(voltage), digits = 1)
    (model == ss) && (plot!(xlabel = "Soma membrane potential (mV)"))
    annotate!([(-40, s / 1.5, text(string(rate) * "Hz", :bottom, :center, :black))])
    annotate!([(-90, s / 1.1, text(string(prox), :bottom, :center, :black, 15))])
    push!(models_hist, p)
end
for x = 1:3
    plot!(models_hist[x], xticks = :none)
end
plot!(models_hist[end], ylabel = "Fluctuation size")
p = plot(
    models_hist...,
    yaxis = false,
    tickfontsize = 13,
    rotation = -45,
    guidefontsize = 18,
    layout = (4, 1),
)

models_hist = []
simtime = 100000
_rate = 1000

for dist = 0:2:6
    model = (150, 350)
    eiratio = [0.0, 0.5, 1.0]
    inputs = TN.make_spikes(
        simtime,
        β,
        baserate = 10.0,
        rate = _rate,
        soma = false,
        dends = true,
        gsyn = 1,
        eiratio = eiratio,
    )

    voltage = run_tripod(model, inputs, simtime, synapses = [0.0, 1.0, dist, 0.0, 1, dist])

    _p = histogram(
        voltage[1, :],
        normalize = true,
        bins = range(-90, stop = -30, length = 180),
        color = :black,
        label = "",
        title = "",
        titlefontsize = 15,
    )
    s = maximum(map(x -> !isnan(x) ? x : 0, _p.series_list[1][:y]))
    kde = globalKDE(3, voltage[1, :])
    plot!(-90:-40, kde, label = "", lw = 2)
    rate = round(TN.get_spike_rate(voltage), digits = 1)
    (model == ss) && (plot!(xlabel = "Soma membrane potential (mV)"))
    annotate!([(-40, s / 1.5, text(string(rate) * "Hz", :bottom, :center, :black))])
    annotate!([(-90, s / 1.1, text(string(dist), :bottom, :center, :black, 15))])
    push!(models_hist, _p)
end
for x = 1:3
    plot!(models_hist[x], xticks = :none)
end
plot!(models_hist[end], ylabel = "Fluctuation size")
q = plot(
    models_hist...,
    yaxis = false,
    tickfontsize = 13,
    rotation = -45,
    guidefontsize = 18,
    layout = (4, 1),
)
z = plot(p, q, xrotation = -45, yticks = :none)
savefig(z, joinpath(@__DIR__, "proximal_distal.pdf"))
z
##

models_hist = []
β = 0.1
simtime = 100000
_rate = 1000

for current = 0:200:600
    model = (350, 150)
    inputs = TN.make_spikes(
        simtime,
        β,
        baserate = 10.0,
        rate = _rate,
        soma = false,
        dends = true,
        gsyn = 1,
        eiratio = [0.0, 1.0, 0.5],
    )

    voltage = run_tripod(model, inputs, simtime, ext_currents = [current, 0.0, 0.0])

    _p = histogram(
        voltage[1, :],
        normalize = true,
        bins = range(-90, stop = -30, length = 180),
        color = :black,
        label = "",
        title = "",
        titlefontsize = 15,
    )
    s = maximum(map(x -> !isnan(x) ? x : 0, _p.series_list[1][:y]))
    kde = globalKDE(3, voltage[1, :])
    plot!(-90:-40, kde, label = "", lw = 2)
    rate = round(TN.get_spike_rate(voltage), digits = 1)
    (model == ss) && (plot!(xlabel = "Soma membrane potential (mV)"))
    annotate!([(-33, s / 1.5, text(string(rate) * "Hz", :bottom, :center, :black))])
    annotate!([(-90, s / 1.2, text(string(current), :bottom, :center, :black, 15))])
    push!(models_hist, _p)
end


for x = 1:3
    plot!(models_hist[x], xticks = :none)
end
plot!(models_hist[end], ylabel = "Fluctuation size")
p = plot(
    models_hist...,
    yaxis = false,
    tickfontsize = 13,
    rotation = -45,
    guidefontsize = 18,
    layout = (4, 1),
)


##

models_hist = []
β = 0.1
simtime = 100000
_rate = 1000

for inh = 0.5:1.0:4.0
    model = (350, 150)
    inputs = TN.make_spikes(
        simtime,
        β,
        baserate = 10.0,
        rate = _rate,
        soma = false,
        dends = true,
        gsyn = 1,
        eiratio = [0.0, 1.0, 0.5],
    )

    voltage = run_tripod(model, inputs, simtime, synapses = [0.0, 1, 1, 0.0, 1.0, 1.0])

    _p = histogram(
        voltage[1, :],
        normalize = true,
        bins = range(-90, stop = -30, length = 180),
        color = :black,
        label = "",
        title = "",
        titlefontsize = 15,
    )
    s = maximum(map(x -> !isnan(x) ? x : 0, p.series_list[1][:y]))
    kde = globalKDE(3, voltage[1, :])
    plot!(-90:-40, kde, label = "", lw = 2)
    rate = round(TN.get_spike_rate(voltage), digits = 1)
    (model == ss) && (plot!(xlabel = "Soma membrane potential (mV)"))
    annotate!([(-40, s / 1.5, text(string(rate) * "Hz", :bottom, :center, :black))])
    annotate!([(-90, s / 1.1, text(string(inh), :bottom, :center, :black, 15))])
    plot!(ylims = (0, s))
    push!(models_hist, _p)
end
for x = 1:3
    plot!(models_hist[x], xticks = :none)
end
plot!(models_hist[end], ylabel = "Fluctuation size")
q = plot(
    models_hist...,
    yaxis = false,
    tickfontsize = 13,
    rotation = -45,
    guidefontsize = 18,
    layout = (4, 1),
)
z = plot(p, q, xrotation = -45, yticks = :none)
savefig(z, joinpath(@__DIR__, "currents.pdf"))
z
