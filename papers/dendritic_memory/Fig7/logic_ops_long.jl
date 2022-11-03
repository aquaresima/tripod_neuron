"""
In this file we will explore the computational property of the Tripod in distinguishing between different input conditions.
Namely:

        A   B   AND OR  XOR IM
_ _     o   o   o   o   o   x
A _     x   o   o   o   x   x
_ B     o   x   o   x   x   o
A B     x   x   x   x   o   x

The neuron activity (membrane, adaptation current) is divided in blocks of 50ms and each block is scored against the
operators. For each operator there is a trained classifier that evaluate the neuron is
the accepted condition.
"""
##

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using StatsBase, Statistics, HDF5 
using MLDataUtils, MLJLinearModels, RollingFunctions
using ProgressBars, Logging

include("stimuli.jl")

function _trainclassifier(;OP::Symbol, X::Matrix{Float64},inputs::Vector{Symbol}, λ=0.5::Float64)
    y = ones(size(inputs))*-1
    for i in eachindex(inputs)
        (inputs[i] ∈ getfield(operators,OP)) && (y[i] = 1)
    end
    train_std = StatsBase.fit(ZScoreTransform, X, dims=2)
    StatsBase.transform!(train_std,X)
    intercept = true
    # deploy MultinomialRegression from MLJLinearModels, λ being the strenght of the reguliser
    lr = LogisticRegression(λ; fit_intercept=intercept)
    # Fit the model
    θ  = MLJLinearModels.fit(lr, X', y)
    return (θ, train_std)
end

function _testclassifier(;OP::Symbol, X::Matrix{Float64},inputs::Vector{Symbol}, lr)
    y = ones(size(inputs))*-1
    for i in eachindex(inputs)
        (inputs[i] ∈ getfield(operators,OP)) && (y[i] = 1)
    end
    θ, train_std = lr
    StatsBase.transform!(train_std,X)
    preds = MLJLinearModels.apply_X(X',θ)
    # #and evaluate the model over the labels
    score = mean(sign.(preds) .== y)
    return score, [inputs, sign.(preds), y]
end

function train_classifiers(X::Matrix{Float64}, inputs::Vector{Symbol})
    classifiers = []
    _ops = []
    for op in keys(operators)
        lr = _trainclassifier(OP=op, X=X,inputs=inputs; λ=0.5::Float64)
        push!(classifiers, lr)
        push!(_ops, op)
    end
    return (;zip(_ops,classifiers)...)
end

function test_classifiers(X::Matrix{Float64}, inputs::Vector{Symbol}, classifiers::NamedTuple)
    scores = []
    confusion = []
    for op in keys(operators)
        score, data = _testclassifier(OP=op, X=X,inputs=inputs,lr=getfield(classifiers,op))
        push!(scores, score)
        push!(confusion, data)
    end
    return scores, confusion
end

function test_model(model, n_symbols, symbol_time)
    @info "Testing model: ", model
    data, inputs = generate_sequence(model, n_symbols, symbol_time)
    features = feats_from_data(data, symbol_time)
    classifiers = train_classifiers(features, inputs)

    data, inputs = generate_sequence(model, n_symbols, symbol_time)
    features = feats_from_data(data, symbol_time)
    scores, confusion = test_classifiers(features, inputs, classifiers)
    return scores, confusion
end


## 
symbol_time = 50
n_symbols = 5000
models = TN.models

# Run simulation
scores = Dict{}()
confusions = Dict{}()
Threads.@threads for m in eachindex(models)
    model = models[m]
    label=TN.labels[m]
    score, confusion = test_model(model, n_symbols , symbol_time )
    print( label,":\n", round.(score, digits=3),"\n")
    push!(scores,label => score )
    push!(confusions, label =>confusion )
end
#

TN.labels
file =datadir("dendritic_memory","logic_ops_data.jld2")
data = @strdict confusions scores symbol_time n_symbols conditions operators
safesave(file, data)
