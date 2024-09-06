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
    plot(
        pp,
        aph,
        layout = (1, 2),
        right_margin = 5Plots.mm,
        bottom_margin = 0Plots.mm,
        top_margin = 0Plots.mm,
    )


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
    plot(NMDAplots[1:4]..., layout = layout, legend = false),
    plot(AMPAplots[1:4]..., layout = layout, legend = false),
    legend = false,
    tickfontsize = 11,
    size = one_column_size .* (1.2, 1.8),
)

@info TN.balance_kie_rate.nmda[10, 10]
p

path = mkpath(plotsdir("up_down", "Figures"))
savefig(p, joinpath(path, "Fig5.pdf"))
