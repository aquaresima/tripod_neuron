##
# Load Tripod
cd(@__DIR__)
using Main
module TripodNeuron
using Logging, Printf, Distributions
# include("TripodNeuron/base/structs/dendrites.jl")
include("../../TripodNeuron/base.jl")
postspike = PostSpike(10, round(Int, 1 / dt), 1 / 50)
AdEx = get_bonoclopath_params()
end

TN = TripodNeuron
using Plots
using LinearAlgebra
using StatsBase
##

words = ["that"]
target_dialects = [1]
target_gender = "m"
using DataFramesMeta
in_words(df_words) = !isempty(intersect(Set(df_words), Set(words)))
in_dialect(df_dialect) = df_dialect ∈ target_dialects
in_gender(df_gender) = occursin(df_gender, target_gender)
hello = @where(train, in_dialect.(:dialect), in_gender.(:gender), in_words.(:words))

spikes, durations = SpikeTimit.get_word_spikes(hello, words[1])


##

n_spikes = round(sum(length.(spikes[1])) / durations[1], digits = 2)
ax = SpikeTimit.raster_plot(spikes[1]; title = "rate: $n_spikes Hz").plt

SpikeTimit.resample_spikes!(spikes, n_feat = 8)

n_spikes = round(sum(length.(spikes[1])) / durations[1], digits = 2)
bx = SpikeTimit.raster_plot(spikes[1], title = "rate: $n_spikes Hz").plt
plot(ax, bx)

SpikeTimit.transform_into_bursts!(spikes)

n_spikes = round(sum(length.(spikes[1])) / durations[1], digits = 2)
cx = SpikeTimit.raster_plot(spikes[1], title = "rate: $n_spikes Hz").plt

##
p = plot(ax, bx, cx, layout = (3, 1))
savefig(p, joinpath(@__DIR__, "results", "raster_plots.pdf"))

##
# n_spikes = sum.(findall.(x-> x>0., spikes[1]))
n_spikes = (length.(spikes[1]))
interval = durations[1]
p = bar(n_spikes, title = "All bins:" * words[1], yticks = :none)
plots = []
for x = 1:4
    push!(plots, bar(n_spikes[x:4:end], title = "Bin: " * string(x)))
end

plot!(plots[1], ylabel = "Input rate (Hz)")

l = @layout [
    c{0.3h}
    grid(2, 2)
]
#
p = plot(
    p,
    plots...,
    layout = l,
    labels = "",
    yticks = (0:10:50, round.(collect(0:10:50) ./ interval, digits = 1)),
)
# savefig(p,  joinpath(@__DIR__,"results","input_rate_$(words[1]).pdf"))

##

words = [
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
target_dialects = [1, 2]
target_gender = "m"
using DataFramesMeta
in_words(df_words) = !isempty(intersect(Set(df_words), Set(words)))
in_dialect(df_dialect) = df_dialect ∈ target_dialects
in_gender(df_gender) = occursin(df_gender, target_gender)
hello = @where(train, in_dialect.(:dialect), in_gender.(:gender), in_words.(:words))
# @show hello
##
#===================================
# Import all the features for the selected words (-> single dialect, male)
and measure the statistical properties.
1. Average rate
2. Average CV ISI
3. PCA
===================================#

in_words.(hello.words)
labels = Vector{String}()
features = Dict(:rate => Vector{Vector{Int}}(), :cv => Vector{Vector{Float32}}())
durations = Vector{Float32}()
hello = @where(train, in_dialect.(:dialect), in_gender.(:gender), in_words.(:words))
for word in words
    if !isempty(SpikeTimit.find_word(; df = hello, words = word))
        global hello
        spikes, intervals = SpikeTimit.get_word_spikes(hello, word)
        SpikeTimit.resample_spikes!(spikes, n_feat = 8)
        SpikeTimit.transform_into_bursts!(spikes)
        for (spike, interval) in zip(spikes, intervals)
            # feats = sum.(findall.(x-> x>0., spike))
            _rate = length.(spike)
            _cv = TN.CV.(diff.(spike))
            push!(labels, word)
            push!(features[:rate], _rate)
            push!(features[:cv], _cv)
            push!(durations, interval)
        end
    end
end
features
labels
## Features and classes

Xtr = vcat(hcat(features[:rate]...), hcat(features[:cv]...))
# Xtr = float.(hcat(features...))
colors = cgrad(:lightrainbow)[0:1/17:17/17]
color_labels = []
correct_class = []
for word in labels
    push!(color_labels, colors[findfirst(x -> x == word, collect(Set(labels)))])
    push!(correct_class, findfirst(x -> x == word, collect(Set(labels))))
end

average_rate = zeros(17, 80)
average_cv = zeros(17, 80)

for class = 1:17
    items = findall(x -> x == class, correct_class)
    average_rate[class, :] .= mean(Xtr[1:80, class], dims = 2)[:, 1]
    average_cv[class, :] .= mean(Xtr[81:end, class], dims = 2)[:, 1]
end

p = plot(
    heatmap(average_rate, title = "Rate", xticks = :none),
    heatmap(average_cv, title = "CV ISI", xlabel = "Input channels"),
    layout = (2, 1),
    yticks = (1:17, words),
    tickfontsize = 8,
)
# savefig(p,  joinpath(@__DIR__,"results","features_space.pdf"))

##
# n_spikes = sum.(findall.(x-> x>0., spikes[1]))
n_spikes = mean(Xtr, dims = 2)[1:80, 1]

p = bar(n_spikes, title = "All bins: all_words", yticks = :none)
plots = []
for x = 1:4
    push!(plots, bar(n_spikes[x:4:end], title = "Bin: " * string(x)))
end

plot!(plots[1], ylabel = "Input rate (Hz)")

l = @layout [
    c{0.3h}
    grid(2, 2)
]
#
p = plot(p, plots..., layout = l, labels = "")
savefig(p, joinpath(@__DIR__, "results", "input_rate_all_words.pdf"))

## GMM
using MultivariateStats
import MultivariateStats: transform
using ColorSchemes

train_std = StatsBase.fit(ZScoreTransform, Xtr, dims = 2)
StatsBase.transform!(train_std, Xtr)
# train a PCA model
M = fit(PCA, Xtr)
# apply PCA model to testing set
Yte = transform(M, Xtr)

μs = zeros(10, 17)

for n = 1:17
    μs[:, n] = mean(Yte[1:10, findall(x -> x == n, correct_class)], dims = 2)
end

gr()
# plotlyjs()
p = scatter(
    Yte[1, :],
    Yte[2, :],
    Yte[3, :],
    c = color_labels,
    lc = :match,
    ms = 4,
    label = "",
)
p = scatter!(μs[1, :], μs[2, :], μs[3, :], c = colors, ms = 10, labels = "Baricenters")
savefig(p, joinpath(@__DIR__, "results", "baricenters.pdf"))

##
using GaussianMixtures

X = copy(Yte[1:10, :]')
gmm = GMM(17, X; method = :kmeans, kind = :full, nInit = 50, nIter = 10, nFinal = 50)
em!(gmm, X)

# scatter!(gmm.μ[:,1], gmm.μ[:,2], gmm.μ[:,3], c=colors, label="", ms=10)
# scatter!(gmm.μ[:,1], gmm.μ[:,2], label="", ms=10)
posterior, lgl = gmmposterior(gmm, X)

colors = cgrad(:lightrainbow)[0:1/17:17/17]
train_class = []
for x in argmax(posterior, dims = 2)
    push!(train_class, x[2])
end

p = plot(
    scatter(Yte[1, :], Yte[2, :], c = colors[train_class], label = ""),
    scatter(Yte[1, :], Yte[2, :], c = color_labels, label = ""),
    markerstrokealpha = 0.0,
    title = "Gaussian Mixtures model",
)
# savefig(p,  joinpath(@__DIR__,"results","GMM clusters.pdf")")


##
function compute_error(target, predict)
    return sum(predict .== target) / length(target)
end


scores = zeros(17, 17)
for l = 1:17
    items = findall(x -> x == l, train_class)
    y = correct_class[items]
    for predict = 1:17
        scores[l, predict] = compute_error(y, predict)
    end
end
argmax(scores, dims = 2)
scatter(maximum(scores, dims = 2))
gmm_score = round(mean(maximum(scores, dims = 2)), digits = 2)
p = heatmap(
    scores,
    title = "Mean GMM score: $gmm_score",
    yticks = (1:17, words),
    xlabel = "Predicted class",
    ylabel = "Matching class",
)
# savefig(p,  joinpath(@__DIR__,"results","GMM scores.pdf"))

##
include("multinomial_regression.jl")
function MultinomialLogisticRegression(
    X::Array{Float32,2},
    labels::Vector{String};
    λ = 0.5::Float32,
    test_ratio = 0.7,
    normalization::Bool = true,
)
    n_classes = length(Set(labels))
    y, _ = string_to_int(labels)
    n_features = size(X, 1)

    # Shuffle the data
    inds = randperm(size(y, 1))
    X = X[:, inds]
    y = y[inds]

    train, test = make_set_index(length(y), test_ratio)

    n_classes_in_train = length(Set(y[train]))
    n_classes_in_test = length(Set(y[test]))
    while !(n_classes_in_train == n_classes && n_classes_in_test == n_classes)
        println("Resampling train and test set, classes were not distributed correctly.")
        train, test = make_set_index(length(y), test_ratio)
        n_classes_in_train = length(Set(y[train]))
        n_classes_in_test = length(Set(y[test]))
    end

    # normalization
    ## Verify the train dataset is normalized
    if normalization
        train_std = StatsBase.fit(ZScoreTransform, X[:, train], dims = 2)
        StatsBase.transform!(train_std, X)
    end

    intercept = false

    # deploy MultinomialRegression from MLJLinearModels, λ being the strenght of the reguliser
    mnr = MultinomialRegression(λ; fit_intercept = intercept)
    # Fit the model
    X_train = copy(X[:, train]')

    θ = MLJLinearModels.fit(mnr, X_train, y[train])
    # # The model parameters are organized such we can apply X⋅θ, the following is only to clarify
    # Get the predictions X⋅θ and map each vector to its maximal element
    # return θ, X
    preds = softmax(MLJLinearModels.apply_X(X[:, test]', θ, n_classes))
    targets = map(x -> argmax(x), eachrow(preds))
    #and evaluate the model over the labels
    y_test = y[test]
    raw_scores = mean(targets .== y_test)
    random_scores = mean(targets .== y_test[randperm(size(y_test, 1))])
    corrected_scores = (raw_scores - random_scores) ./ (1 .- random_scores)

    params = reshape(θ, n_features + Int(intercept), n_classes)

    return corrected_scores, random_scores, params
end

scores = zeros(4, 20)
for x = 1:20
    # score = -1
    for (z, data) in zip(1:4, [Xtr, Xtr[1:80, :], Xtr[81:160, :], Yte])
        score = 0.0
        while score <= 0.1
            score, randomscore, A = MultinomialLogisticRegression(data, labels)
            scores[z, x] = score
        end
    end
end
p = bar(
    mean(scores, dims = 2),
    yerror = std(scores, dims = 2),
    ylabel = "LogReg classification",
    xticks = (1:4, ["Rate + CV", "Rate", "CV", "PCA"]),
    label = false,
)
# savefig(p,  joinpath(@__DIR__,"results","logreg_input.pdf"))


# histogram(A[:,1:3], labels="", title="score: "* string(score), bins = -2:0.3:2)
