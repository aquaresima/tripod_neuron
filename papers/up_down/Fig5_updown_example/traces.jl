using DrWatson
using Plots, Revise, Random
@quickactivate "Tripod"
using TripodNeuron
using Statistics
include(projectdir("scripts", "up_down", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_plots.jl"));

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
        if β == 0
            plot!(pp, title = "NMDA + AMPA")
        end
    else
        c = :darkred
        if β == 0
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
    )
    aph = APH(voltage, "", filter = false)
    plot!(
        aph,
        yticks = false,
        yaxis = false,
        frame = :axes,
        xlabel = "Soma (mV)",
        xrotation = -45,
        title = "β = $β",
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

TN.μmem
rate = 11
@info TN.νs[rate]
model = (ds = (400, 400), nmda = true)
NMDAplots = []
for β = 0:80:240
    p = get_trace(model, β, rate, do_spikes = true)
    push!(NMDAplots, p)
end
model = (ds = (400, 400), nmda = false)
AMPAplots = []
for β = 0:80:240
    p = get_trace(model, β, rate, do_spikes = true)
    push!(AMPAplots, p)
end

# TN.ampa_equivalent.Esyn_dend.AMPA.g0
for (p, q) in zip(NMDAplots[1:end-1], AMPAplots[1:end-1])
    plot!(p, frame = :none, xlabel = "")
    plot!(q, frame = :none, xlabel = "")
end

layout = @layout[
    a{0.25h}
    a{0.25h}
    a{0.25h}
    a{0.25h}
]
p = plot(
    plot(NMDAplots[1:4]..., layout = layout, legend = false, size = (650, 800)),
    plot(AMPAplots[1:4]..., layout = layout, legend = false, size = (650, 800)),
    legend = false,
    tickfontsize = 11,
)

@info TN.balance_kie_rate.nmda[10, 10]
p

##
all_spikes = []
for β = 1:1:240
    for _ in [1, 2]
        spikes = TN.make_spikes(1_000, β, rate = 50, seed = 29)[3, :]
        push!(all_spikes, rand(-20:20) .+ (findall(spikes .> 0)) .* TN.dt)
    end
end

all_spikes

raster_plot(reverse(Spiketimes(all_spikes)));
q2 = plot!(
    frame = :axes,
    size = (800, 800),
    yticks = (0:160:480, reverse(collect(0:80:240))),
    xaxis = false,
    ylabel = "Fluctuation size (β)",
    title = "input spike train",
    xticks = (0:250:1000, 0:0.25:1),
    xlabel = "Time (s)",
)


##
r = []
for β = 10:60:250
    push!(r, TN.make_rates(1_000, β, rate = 1000, seed = 29))
end
q3 = plot(
    r,
    c = collect(palette(:roma, 5))',
    lw = 3,
    labels = collect(0:60:240)',
    legendtitle = "β",
    ylabel = "Input rate (kHz)",
    xticks = :none,
    yticks = ((500, 1000, 2000), (0.5, 1, 2)),
    ylims = (500, 2500),
    xlabel = "",
)

# q3 = 
# plot!(ylabel="Input β", xlabel="", xticks=:none )
##
data = load(datadir("up_down", "robustness", "cv_inputs.bson"))[:data]

q1 = plot(
    data.βs,
    mean(data.cvs[:, :]', dims = 2),
    ribbon = std(data.cvs[:, :]', dims = 2),
    label = "",
    xlabel = "Fluctuation size (β)",
    ylabel = "ISI CV",
)
layout = @layout [
    a{0.3h}
    b{0.4h}
    b{0.3h}
]
q = plot(q1, q3, q2, layout = layout)

##
layout = @layout[a{0.3w} b{0.7w}]
z = plot(q, p, layout = layout, size = (1200, 1000))
path = mkpath(plotsdir("up_down", "up_down", "robustness"))
savefig(z, joinpath(path, "input_specifics.pdf"))
z


p
path = mkpath(plotsdir("up_down", "Figures", "Fig5.pdf"))




##For poster
r = []
for β = 10:60:250
    push!(r, TN.make_rates(1_000, β, rate = 1000, seed = 29))
end
q3 = plot(
    r,
    c = collect(palette(:roma, 5))',
    lw = 3,
    labels = collect(0:60:240)',
    legendtitle = "β",
    ylabel = "Input rate (kHz)",
    xlabel = "Time (s)",
    yticks = ((500, 1000, 2000), (0.5, 1, 2)),
    ylims = (500, 2500),
    xlims = (0, 1500),
    xticks = (0:250:1000, 0:0.25:1),
    frame = :axes,
)

# q3 = 
# plot!(ylabel="Input β", xlabel="", xticks=:none )
data = load(datadir("up_down", "robustness", "cv_inputs.bson"))[:data]

q1 = plot(
    data.βs,
    mean(data.cvs[:, :]', dims = 2),
    ribbon = std(data.cvs[:, :]', dims = 2),
    label = "",
    xlabel = "Fluctuation size (β)",
    ylabel = "ISI CV",
)
layout = @layout [
    a{0.3h}
    b{0.4h}
    b{0.3h}
]
q = plot(q1, q3, layout = layout)

path = mkpath(plotsdir("up_down", "up_down", "robustness"))
savefig(q, joinpath(path, "input_specifics_poster.pdf"))
