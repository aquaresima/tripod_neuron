include("learning_functions.jl")
include("get_datasets.jl")
##
using HDF5
using JLD
using StatsBase
using EllipsisNotation
using DataFramesMeta
using ColorSchemes
using Plots

##

## Run with file for different words as argument
IN_NEURONS = 80
@show ARGS
if length(ARGS) > 0
    my_arg = parse(Int, ARGS[1])
    WORD = my_arg:my_arg
else
    WORD = 1:1
end


## Set the input
common_words = [
    "that",
    "had",
    "she",
    "me",
    "your",
    "all",
    "like",
    "don't",
    "year",
    "water",
    "dark",
    "rag",
    "oily",
    "wash",
    "ask",
    "carry",
    "suit",
]
words = common_words[WORD]
target_dialects = [1]
target_gender = "m"

println("Running for word: $words")
in_words(df_words) = !isempty(intersect(Set(df_words), Set(words)))
in_dialect(df_dialect) = df_dialect ∈ target_dialects
in_gender(df_gender) = occursin(df_gender, target_gender)

df = @where(train, in_dialect.(:dialect), in_gender.(:gender), in_words.(:words))


# spikes, durations = SpikeTimit.get_word_spikes(df, words)
# SpikeTimit.resample_spikes!(spikes, n_feat=8)
# SpikeTimit.transform_into_bursts!(spikes)
durations, spikes, labels = SpikeTimit.select_inputs(df, words, samples = 10, n_feat = 8)

## Run the vstdp learning rule
EFF = 3.0
s = TN.AdEx_Soma()
all_ft, all_neurons, all_words, all_phones =
    SpikeTimit.mix_inputs(durations, spikes, labels, repetitions = 5, silence_time = 0.0)
all_ft_dt = SpikeTimit.all_ft_to_dt(all_ft, TN.dt)
SIMTIME = round(Int, ceil(all_ft[end]) * 1000)

AMPAsynapses = EFF .* ones(IN_NEURONS)
NMDAsynapses = EFF .* ones(IN_NEURONS)


ampa_c, voltage_c, _, rate_c = simulate_clopath(
    s,
    SIMTIME,
    NMDAsynapses,
    AMPAsynapses,
    all_ft_dt,
    all_neurons,
    max_efficacy = EFF * 2,
    inh_τ = 1000.0,
)
## Run the triplet learning rule

s = TN.AdEx_Soma()
all_ft, all_neurons, all_words, all_phones =
    SpikeTimit.mix_inputs(durations, spikes, labels, repetitions = 5, silence_time = 0.0)
all_ft_dt = SpikeTimit.all_ft_to_dt(all_ft, TN.dt)
SIMTIME = round(Int, ceil(all_ft[end]) * 1000)
AMPAsynapses = EFF .* ones(IN_NEURONS)
NMDAsynapses = EFF .* ones(IN_NEURONS)

ampa_p, voltage_p, _, rate_p = simulate_pfister(
    s,
    SIMTIME,
    NMDAsynapses,
    AMPAsynapses,
    all_ft_dt,
    all_neurons,
    max_efficacy = EFF * 2,
    inh_τ = 1000.0,
)

## Plot results
tt = length(voltage_c)
rc = round(TN.get_spike_rate(voltage_c), digits = 2)
rp = round(TN.get_spike_rate(voltage_p), digits = 2)
rc = round(TN.CV(voltage_c), digits = 2)
rp = round(TN.CV(voltage_p), digits = 2)


dt = round(tt / 10000, digits = 2)
Plots.default(background_legend = :white)
p = plot(
    plot(voltage_p),
    plot(voltage_c, xlabel = "Time (s)", ylabel = "Membrane\npotential (mV)"),
    title = ["Triplet ($rp (Hz))" "vSTDP ($rc (Hz))"],
    layout = (2, 1),
    labels = "",
    xticks = (0:tt/5:tt, string.(0:dt/5:dt)),
)
annotate!(p, 2000, 20, Text("word: $(words)"))

q = plot(
    plot(ampa_p[:, 1:20:end]', labels = ""),
    plot(
        ampa_c[:, 1:20:end]',
        labels = "",
        xlabel = "Time (s)",
        ylabel = "Synaptic \efficacy (a.u)",
    ),
    title = ["Triplet ($rp (Hz))" "vSTDP ($rc (Hz))"],
    layout = (2, 1),
    labels = "",
    xticks = (0:tt/5/20:tt/20, string.(0:dt/5:dt)),
)
annotate!(q, 2000, 20, Text("word: $(words)"))
savefig(p, joinpath(@__DIR__, "results", "voltage_$(words).pdf"))
savefig(q, joinpath(@__DIR__, "results", "synapses_$(words).pdf"))

# p
# q
##
sum(voltage_c .> 0)
sum(voltage_p .> 0)
plot([voltage_c, voltage_p], xlims = (1, 1e4))
# plot
## Analysis
function plot_synapses_rate(ampas, labels, spikes, durations)

    rate = zeros(80)
    for (spike, dur) in zip(spikes, durations)
        rate .+= (length.(spike) ./ dur)
    end

    plots = []

    for x = 1:4
        rate .= rate ./ maximum(rate)
        p = bar(
            rate[x:4:80],
            alpha = 0.3,
            ylabel = "Input rate",
            ylims = (0, 1),
            yticks = :none,
            label = "",
        )
        p = plot!(twinx(), yticks = :none, ylabel = "Syn. weights")
        # p = scatter!(twinx(), mean_ampa[x:4:80], c=c[sorting[x:4:80]], ylims=(-0.2,max_ampa*1.3), yticks=:none, ylabel="Syn. weights", ms=7, alpha=1)
        for (a, c, l) in zip(ampas, palette(:tab10), labels)
            _l = x == 1 ? l : ""
            _a = mean(a, dims = 2)[:, 1]
            max_ampa = maximum(_a)
            p = scatter!(
                _a[x:4:80] ./ max_ampa,
                c = c,
                ylims = (-0.2, 1.3),
                markerstrokealpha = 0.0,
                yticks = :none,
                ylabel = "Syn. weights",
                ms = 7,
                alpha = 1,
                label = _l,
            )
        end
        plot!(
            p,
            title = "Channel $x",
            legendfontsize = 11,
            legend = :bottomleft,
            background_legend = :transparent,
        )
        # plot()
        # if x == 1
        # 	plot!(p, [[],[]], labels=labels, c=palette(:tab10)[1:2]', legendfontsize=9, legend=:bottomleft)
        # end
        push!(plots, p)
    end
    l = @layout [
        a{0.45w} x{0.1w} b{0.45w}
        c{0.45w} y{0.1w} d{0.45w}
    ]
    insert!(plots, 4, plot(frame = :none, xlabel = "\nInput frequency"))
    insert!(
        plots,
        2,
        plot(frame = :none, title = "Input rates vs synaptic weights \nword: $(words)\n"),
    )
    return plot(plots..., layout = l, guidefontsize = 13)
end

p = plot_synapses_rate([ampa_p, ampa_c], ["triplet", "vstdp"], spikes, durations)
savefig(p, joinpath(@__DIR__, "results", "synapse_rate_$(words).pdf"))

p
# plot(plot_synapses_rate(ampa_p, spikes, durations),plot_synapses_rate(ampa_c, spikes, durations))
##
function compare_synapses(ampa, title, spikes, durations)

    rate = zeros(80)
    for (spike, dur) in zip(spikes, durations)
        rate .+= (length.(spike) ./ dur)
    end


    c = cgrad(:roma)[0:1/80:1]
    mean_ampa = mean(ampa, dims = 2)[:, 1]
    sorting = sort(1:80, by = y -> sort(1:80, by = x -> mean_ampa[x], rev = true)[y])
    # sorting = sort(1:80, by=y->sort(1:80,by=x->mean_ampa[x], rev=true)[y])
    #
    l = @layout [
        a{0.3h}
        b{0.3h}
        c{0.3h}
    ]
    return plot(
        scatter(
            mean(ampa, dims = 2),
            StatsBase.std(ampa, dims = 2),
            label = "",
            c = c[sorting],
            ms = 7,
            ylabel = "Std Dev \nweights",
            title = title,
        ),
        scatter(
            mean(ampa, dims = 2),
            ampa[:, end],
            label = "",
            c = c[sorting],
            ms = 7,
            ylabel = "Final \n weight",
            ylims = (-0.2, maximum(ampa[:, end]) + 0.2),
        ),
        scatter(
            mean(ampa, dims = 2),
            rate,
            ylabel = "Channel \n rate (Hz)",
            xlabel = "Average synaptic weight",
            c = c[sorting],
            label = "",
            ms = 7,
        ),
        layout = l,
        guidefontsize = 14,
        ms = 5,
        markerstrokealpha = 0.0,
    )
end

p = plot(
    compare_synapses(ampa_p, "triplet", spikes, durations),
    plot(compare_synapses(ampa_c, "vSTDP", spikes, durations), ylabel = ""),
)
# savefig(p, joinpath(@__DIR__, "results", "synapse_statistics_*$(words).pdf"))
p
##

function compare_learning(ampa_p, ampa_c, spikes, durations)

    rate = zeros(80)
    for (spike, dur) in zip(spikes, durations)
        rate .+= (length.(spike) ./ dur)
    end
    c = cgrad(:roma)[0:1/80:1]
    sorting = sort(1:80, by = y -> sort(1:80, by = x -> rate[x], rev = true)[y])

    ampa = cat(ampa_p, ampa_c, dims = 3)
    ampa = mean(ampa, dims = 2)[:, 1, :]
    return scatter(
        ampa[:, 1],
        ampa[:, 2],
        c = c[sorting],
        xlabel = "Triplet",
        ylabel = "vSTDP",
        ms = 12,
        labels = "rate (Hz)",
        legend = :bottomright,
    )
end

p = compare_learning(ampa_p, ampa_c, spikes, durations)
# savefig(p, joinpath(@__DIR__, "results", "triplet_vs_vSTDP$(words).pdf"))
p
