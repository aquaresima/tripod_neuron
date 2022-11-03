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


symbol_time = 50
models = TN.models

function get_data_logic(u_r = -70)
    membrane = Dict()
    symbols = Dict()
    features = Dict()

    Threads.@threads for m in eachindex(models)
        model = models[m]
        label=TN.labels[m]
        vw, symbol = generate_sequence(model, conditions, symbol_time,u_r )
        feats = feats_from_data(copy(vw), symbol_time)
        av_feats = zeros(size(feats)...,100)
        for x in 1:100
            vw, symbol = generate_sequence(model, conditions, symbol_time, u_r )
            feats = feats_from_data(copy(vw), symbol_time)
            av_feats[:,:,x] .= feats
        end
        feats = mean(av_feats, dims=3)[:,:,1]


        push!(membrane,label => vw )
        push!(symbols, label => symbol )
        push!(features, label=>feats)
    end
    return features, membrane, symbols
end


features,membrane,symbols= get_data_logic()
LETTERS =  string.(collect('A':'Z')[1:16])
plots = []
for (n, key) in enumerate(TN.labels)
    vw = membrane[key]

    for x in 1:length(conditions)
        title= n ==1 ? ": $(symbols[key][x])" : ""
        ylabel= n ==3  ? "                        Membrane potential (mV)" : ""
        ylabel= x == 1 ? "$key \n \n $ylabel" : ""
        yylabel= x==4 && n ==3 ? "                          Adaptive current (nA)" : ""

        p = plot(vw[1:500,x],ylabel=ylabel, title = title,titlefontsize=18, lw=4, c=:black, alpha=0.5, ylims=(-90,25))

        letter = popfirst!(LETTERS)
        annotate!((10, 35, Plots.text(letter, 18, :black, :bold)))
        # (x >1) && (plot!(p, yticks=:none))
        q = twinx(p)
        q = plot!(q,vw[501:end,x]/1000, c=RED,  ylabel=yylabel, lw=4, ylims=(0,1.5))
        plot!(xticks=(0:250:500,0:25:50), legend=false)
        (n <4) && (plot!(q, xticks=:none, xlabel=""))
        (n==4) && (x==2) && plot!(xlabel="                           Time (ms)")
        push!(plots, q,)
    end
    # push!(plots,plot(frame=:none))
end

# layout= @layout ([grid(4,4), [a{0.1h} b{0.1h} c{0.1h} d{0.1h}]])

layout= @layout ([ grid(4,4){0.90w} [_  _ _ _]])

pp = plot(plots..., size=(1200,1200), layout=layout, 
    right_margin=10Plots.mm, top_margin=5Plots.mm)

savefig(pp, plotsdir("dendritic_memory", "Fig7_appendix.pdf"))

pp