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
fig_path = mkpath(plotsdir("up-down", "Figures"))
##
τ = 50
data = load(datadir("up_down", "bistability", "stimulus_betas_$τ.jld2")) |> dict2ntuple
function plot_cum_prob(data, bs, state, receptor = "NMDA"; color = :black, kwargs...)
    plot()
    @unpack condition, βs = data
    if length(bs) == 2
        for (b, ls) in zip(bs, [:solid :dash :dot])
            states = getfield(condition[b][receptor], state)[1, :] ./ 1000
            m = maximum(states)
            _hist = length.(UncertainData.bin(0.:0.01:m, states, 1:length(states)))
            plot!(0.01:0.01:m, cum_prob(_hist), c = color, lw =3, ls = ls, label = "")
        end
    else
        color = palette(:viridis, length(bs))
        for (n,b) in enumerate(bs)
            states = getfield(condition[b][receptor], state)[1, :] ./ 1000
            m = maximum(states)
            @info m
            _hist = length.(UncertainData.bin(0.:0.01:m, states, 1:length(states)))
            plot!(0.01:0.01:m, cum_prob(_hist), c = color[n], lw =3, label = "")
        end
    end
    plot!(; kwargs...)
end
##
data.βs

# bs = 2:1:11

bs = [2,6]
z1 = plot(
    plot_cum_prob(
        data,
        bs,
        :up,
        "AMPA",
        color = :darkred,
        # legend = false,
        leftmargin = 10Plots.mm,
        # yticks = (0:0.5:1, 0:0.5:1),
    ),
    ylims = (0, 1),
    xlims = (-0.2, 11),
    topmargin = 10Plots.mm,
    yticks = :none,
    title= "UP",
    xticks = :none,
    legend = false,
)

z2 = plot(
    plot_cum_prob(
        data,
        bs,
        :down,
        "AMPA",
        color = :darkred,
        leftmargin = 10Plots.mm,
    ),
    layout = (1, 3),
    ylims = (0, 1),
    xlims = (-0.2, 2),
    bottommargin = 10Plots.mm,
    title= "DOWN",
    topmargin = 3mm,
    xticks = :none,
    yticks = 0:0.5:1,
    # yticks = (0:0.2:1, 0:20:100),
)
beta1, beta2 = [data.βs[1], data.βs[5] ]./1000 .+1
plot!(
    [[], []],
    ls = [:solid :dash],
    # label = ["ISI CV = $beta1" "ISI CV = $beta2"],
    label = [L"ISI_{CV}" *" $beta1" L"ISI_{CV}"*" $beta2"],
    c = :black,
    lw = 1,
)


plot!(z1, bottommargin=0Plots.mm)# legend = false)
plot!(z2, bottommargin=0Plots.mm)# legend=false)

z3 = plot(
    plot_cum_prob(
        data,
        bs,
        :up,
        "NMDA",
        color = :darkblue,
        # legend = false,
        leftmargin = 10Plots.mm,
        # yticks = (0:0.5:1, 0:0.5:1),
    ),
    xticks = (1:2:11),
    ylims = (0, 1),
    xlims = (-0.2, 11),
    topmargin = 10Plots.mm,
    yticks = :none,
)


z4 = plot(
    plot_cum_prob(
        data,
        bs,
        :down,
        "NMDA",
        color = :darkblue,
        legend = false,
        ylabel="           Cumulative probability",
        leftmargin = 10Plots.mm,
    ),
    layout = (1, 3),
    # legend = false,
    xlabel = "Time (s)",
    xticks = (0:1:3),
    ylims = (0, 1),
    xlims = (-0.2, 2),
    yticks = 0:0.5:1,
    bottommargin = 10Plots.mm,
    topmargin = 3mm,
)
plot!(z3, topmargin=0Plots.mm, legend = false)
plot!(z4, topmargin=0Plots.mm, legend=false)
zz1 = plot(z2, z1, z4, z3, layout = (2, 2), size = (400, 500), leftmargin=3Plots.mm, rightmargin=3Plots.mm, frame=:axes, legendfontsize=10)

savefig(zz1, joinpath(fig_path, "Fig7_duration.pdf"))
zz1
# plot!(zz1, subplot=2, rightmargin=10Plots.mm)
# plot!(zz1, subplot=1, leftmargin=10Plots.mm, size=(600,400))
##
β= 100
simtime = 10_000
model = TN.models[3]
_rate = 11
# τ = 50
@unpack νs = TN
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
_rate = 11;
@show νs[_rate];
plots = map(0:3) do n
    τ = 5^n
    Random.seed!(133)
    samples = map(1:50) do x
        cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
        inputs = TN.make_spikes(simtime, β; cond..., τ = τ, seed=x)
        voltage =
            TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
        v = voltage[1, 1:100:end]
        v[v.>AdEx.θ] .= AdEx.θ
        v
    end
    vv = vcat(samples...)
    n=round(Int,2*length(vv)/50)
    p1 = plot(vv[:50:n], label="", c=:black, xticks=:none, bottommargin=0Plots.mm)
end
plot!(plots[3], ylabel="Membrane potential (mV)")
[plot!(p, topmargin=0Plots.mm, bottommargin=0Plots.mm) for p in plots]
plot!(plots[end],xticks=(0:50_0:simtime, 0:2.5:10), xlabel="Time (s)", bottommargin=5Plots.mm)
p1= plot(plots..., layout=(4,1), legend=false, size=(600,400), yticks=(-70:10:-50))

##
function residual_probability(bs, data)
    # histogram(data.condition[2]["NMDA"][:up][1,:]./1000, bins=0:1:50)
    nmda_y = map(bs) do x
        states = data.condition[x]["NMDA"][:up][1,:]./1000 
        m = maximum(states)
        durations = 0:0.01:m
        _hist = length.(UncertainData.bin(durations, states, 1:length(states)))
        x = findfirst(durations .== 1)
        1-cum_prob(_hist)[x]    
    end
    ampa_y = map(bs) do x
        states = data.condition[x]["AMPA"][:up][1,:]./1000 
        m = maximum(states)
        durations = 0:0.01:m
        _hist = length.(UncertainData.bin(durations, states, 1:length(states)))
        x = findfirst(durations .== 1)
        1-cum_prob(_hist)[x]    
    end
    p1 = scatter(data.βs[bs]./1000 .+1, nmda_y, smooth=true, c=:darkblue, msc=:darkblue, ms=5, label="NMDA")
    scatter!(data.βs[bs]./1000 .+1,ampa_y,smooth=true, c=:darkred, ms=5, msc=:darkred, shape=:diamond, label="AMPA")
    plot!(yticks=:none, frame=:box, legend=:top, ylims=(0,1))
    plot!(rightmargin=0Plots.mm)
    return p1
end

τs = [1,5,25,50]
bs = 2:11
plots = map(τs) do τ
    data = load(datadir("up_down", "bistability", "stimulus_betas_$τ.jld2")) |> dict2ntuple
    p = residual_probability(bs, data)
    p
end
plot!(plots[1], ylabel="Prob. interval > 1 s ", yticks=(0:0.25:1))
plot!(plots[2], ylabel="", legend=false)
plot!(plots[3], ylabel="", xlabel="ISI CV", legend=false)
plot!(plots[4], ylabel="", legend=false)
p2 = plot(plots..., layout = (1,4), size=(600,600), xrotation=45, bottommargin=8Plots.mm, leftmargin=0Plots.mm, rightmargin=0Plots.mm, tickfontsize=11)

plot!(p1, rightmargin=20Plots.mm)

plots = map(τs) do τ
    data = load(datadir("up_down", "bistability", "stimulus_betas_$τ.jld2")) |> dict2ntuple
    plot(
        plot_cum_prob(
            data,
            bs,
            :up,
            "NMDA",
            color = :darkblue,
            # legend = false,
            leftmargin = 10Plots.mm,
            # yticks = (0:0.5:1, 0:0.5:1),
        ),
        xticks = (1:2:11),
        ylims = (0, 1),
        xlims = (-0.2, 11),
        topmargin = 5Plots.mm,
        yticks = :none,
    )
    annotate!((0.5,1.1),text("τ:$τ ms"))
end
plot!(plots[3], xlabel="State duration (s)")
plot!(plots[1], ylabel="Cumulative prob.", yticks=(0:0.25:1))
p3 = plot(plots..., layout = (1,4), size=(600,600), xrotation=45, bottommargin=8Plots.mm, leftmargin=0Plots.mm, rightmargin=0Plots.mm, tickfontsize=11)

layout = @layout [a{0.45w} [ b{0.5h, 0.95w} 
                            c{0.5h, 0.95w}] 
                ]
zz2 = plot(p1,p3,p2, layout=layout, size=(1200,600))
savefig(zz2, joinpath(fig_path, "Fig7_tau_duration.pdf"))