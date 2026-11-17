## Load learning and Tripod module

include("get_datasets.jl")
include("learning_functions.jl")

## Set input space and prepare spikes

using EllipsisNotation
using DataFramesMeta
IN_NEURONS = 80
WORD = 1:5

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

in_words(df_words) = !isempty(intersect(Set(df_words), Set(words)))
in_dialect(df_dialect) = df_dialect ∈ target_dialects
in_gender(df_gender) = occursin(df_gender, target_gender)

df = @where(train, in_dialect.(:dialect), in_gender.(:gender), in_words.(:words))

spikes, durations = SpikeTimit.get_word_spikes(df, words)
SpikeTimit.resample_spikes!(spikes, n_feat = 8)
SpikeTimit.transform_into_bursts!(spikes)
original_spikes = deepcopy(spikes)

## Define a frequency response function to

function frequency_response(
    inputs,
    frequencies::Vector{Int64},
    learning::Function;
    samples = 50,
    efficacy = 4.0,
    mix_inputs = false,
    kwargs...,
)
    ## Select the frequencies of the input space
    durations, spikes, labels = inputs
    SpikeTimit.select_frequencies!(spikes, frequencies = frequencies)

    ## Prepare spikes
    all_ft, all_neurons, all_words, all_phones, _ids = SpikeTimit.mix_inputs(
        durations,
        spikes,
        labels,
        repetitions = 1,
        silence_time = 0.0,
        ids = true,
    )
    simtime = round(Int, ceil(sum(durations) * 1000))

    ## Run for N samples
    spiketimes = Vector{Vector{Float32}}()
    AMPA = zeros(size(frequencies)[1], simtime, samples)
    ids = zeros(Int, length(durations), samples)
    Threads.@threads for n = 1:samples
        ## Reset the synapses
        AMPAsynapses = efficacy .* ones(IN_NEURONS)
        NMDAsynapses = efficacy .* ones(IN_NEURONS)
        s = TN.AdEx_Soma()
        if mix_inputs
            all_ft, all_neurons, all_words, all_phones, _ids = SpikeTimit.mix_inputs(
                durations,
                spikes,
                labels,
                repetitions = 1,
                silence_time = 0.0,
                ids = true,
            )
        end
        _ft = deepcopy(SpikeTimit.all_ft_to_dt(all_ft, TN.dt))
        _neurons = deepcopy(all_neurons)
        ampa, voltage, θ, rate = learning(
            s,
            simtime,
            NMDAsynapses,
            AMPAsynapses,
            _ft,
            _neurons,
            max_efficacy = efficacy * 2,
            inh_τ = 1000.0,
        )

        ids[:, n] .= _ids
        AMPA[:, :, n] .= ampa[frequencies, 1:10:end]
        push!(spiketimes, TN.get_spike_times(voltage) .* 0.1)
    end
    return AMPA, spiketimes, ids
end

##
using HDF5
using JLD
using StatsBase
using Random


SAMPLES = 30
EFFICACY = 6.5
ITEMS_PER_WORD = 20

inputs = SpikeTimit.select_inputs(df, words, samples = ITEMS_PER_WORD, n_feat = 8)

function frequency_response_run(samples, efficacy, mix_inputs::Bool)
    for x = 1:4
        frequencies = collect(1+(x-1)*20:x*20)
        _inputs = deepcopy(inputs)

        println("Frequency set: $x, mix_input: $mix_inputs ")
        if mix_inputs
            file = joinpath(@__DIR__, "data", "ids_mixed_data_$x.h5")
        else
            file = joinpath(@__DIR__, "data", "ids_data_$x.h5")
        end

        isfile(file) ? rm(file) : nothing
        fid = h5open(file, "w")
        fid["frequences"] = frequencies
        fid["durations"] = durations

        # run for triplet
        myampa, myspikes, ids = frequency_response(
            _inputs,
            frequencies,
            simulate_pfister,
            samples = samples,
            efficacy = efficacy,
            mix_inputs = mix_inputs,
        )
        M = zeros(length(myspikes), maximum(length.(myspikes)))
        for s in eachindex(myspikes)
            M[s, 1:length(myspikes[s])] .= myspikes[s]
        end
        _vstdp = create_group(fid, "vstdp")
        _vstdp["ampa"] = myampa
        _vstdp["spikes"] = M
        _vstdp["ids"] = ids

        # run for vSTDP
        myampa, myspikes, ids = frequency_response(
            _inputs,
            frequencies,
            simulate_clopath,
            samples = samples,
            efficacy = efficacy,
            mix_inputs = mix_inputs,
        )
        M = zeros(length(myspikes), maximum(length.(myspikes)))
        for s in eachindex(myspikes)
            M[s, 1:length(myspikes[s])] .= myspikes[s]
        end
        _triplet = create_group(fid, "triplet")
        _triplet["ampa"] = myampa
        _triplet["spikes"] = M
        _triplet["ids"] = ids
        close(fid)
    end
end

##
length.(myspikes)
my_ampa_c, myspikes_c =
    frequency_response(inputs, collect(40:60), simulate_clopath, samples = 50)

function synapse_stability(synapses, init = 2_000)
    syn = (@view synapses[:, init:end, :])
    return StatsBase.mean(mean(syn, dims = 3), dims = 2)[:, 1, 1],
    StatsBase.std(mean(syn, dims = 3), dims = 2)[:, 1, 1]
end
μ, σ = synapse_stability(my_ampa)

rate = zeros(80)
for (spike, dur) in zip(spikes, durations)
    rate .+= (length.(spike) ./ dur)

end

using StatsPlots
using StatsBase
ctg = repeat(["Category 1", "Category 2"], inner = 5)
nam = repeat("G" .* string.(1:5), outer = 2)

groupedbar(
    nam,
    rand(5, 2),
    group = ctg,
    xlabel = "Groups",
    ylabel = "Scores",
    title = "Scores by group and category",
    bar_width = 0.67,
    lw = 0,
    framestyle = :box,
)
ctg = repeat(["Rate", "Triplet"], inner = 20)
name = repeat(1:20, outer = 2)
groupedbar(
    name,
    vcat(rate[frequencies] ./ maximum(rate) * maximum(μ), μ),
    group = ctg,
    lw = 0.0,
    yerror = hcat(zeros(20), σ),
    bar_width = 0.67,
)# yerror=[zeros(20) σ])
##
TN.raster_plot(myspikes_c)
## Select input

x = 2
using HDF5
using StatsBase
# for x in 1:4
file = joinpath(@__DIR__, "data", "ids_data_$x.h5")
fid = read(h5open(file, "r"))

frequencies = fid["frequences"]
durations = fid["durations"]
clopath = fid["vstdp"]
pfister = fid["triplet"]

p = plot(mean(pfister["ampa"], dims = 3)[:, :, 1]', labels = "")
q = plot(mean(clopath["ampa"], dims = 3)[:, :, 1]', labels = "")

file = joinpath(@__DIR__, "data", "ids_mixed_data_$x.h5")

frequencies = fid["frequences"]
durations = fid["durations"]
clopath = fid["vstdp"]
pfister = fid["triplet"]

r = plot(mean(pfister["ampa"], dims = 3)[:, :, 1]', labels = "", title = "Triplet")
s = plot(mean(clopath["ampa"], dims = 3)[:, :, 1]', labels = "", title = "vSTDP")
plot(r, s)
plot(p, q, r, s)
##

# matrix2st(M::Matrix{Float32}) = map(x->x[x.>0.], eachrow(M))
# matrix2st(pfister["spikes"])

# pfister["ampa"]


# [:,1,1],	ribbons=StatsBase.std(mean(pfister["ampa"], dims=2)[:,1,:],dims=2), labels="")
# SpikeTimit.raster_plot(matrix2st(pfister["spikes"])).plt

plot()
plot!(
    mean(pfister["ampa"], dims = [3, 2])[:, 1, 1],
    ribbons = StatsBase.std(mean(pfister["ampa"], dims = 2)[:, 1, :], dims = 2),
    labels = "",
)
# end
# plot(mean(clopath["ampa"], dims=[2,1])[1,1,:], labels="")


##
plot(
    plot(
        mean(clopath["ampa"], dims = 3)[1:4:20, :, 1]',
        ribbons = StatsBase.std(clopath["ampa"], dims = 3)[:, :, 1]',
        labels = "",
    ),
    plot(
        mean(pfister["ampa"], dims = 3)[1:4:20, :, 1]',
        ribbons = StatsBase.std(pfister["ampa"], dims = 3)[:, :, 1]',
        labels = "",
    ),
    layout = (2, 1),
)
##
file = joinpath(@__DIR__, "results", "data_$x.h5")
fid = read(h5open(file, "r"))
frequencies = fid["frequences"]
durations = fid["durations"]
clopath = fid["vstdp"]
pfister = fid["triplet"]

matrix2st(M::Matrix{Float32}) = map(x -> x[x.>0.0], eachrow(M))

histogram(pfister["spikes"], dims = 2)
matrix2st(pfister["spikes"])

SpikeTimit.raster_plot(matrix2st(pfister["spikes"])).plt

# durations
