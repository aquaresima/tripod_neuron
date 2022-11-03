"""
This file reproduces Fig 9a
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using StatsPlots, CategoricalArrays, LaTeXStrings, RecipesPipeline, Statistics, JLD2


file =datadir("dendritic_memory","logic_ops_data.jld2")
scores_mem = JLD2.load(file)["scores"]
confusions = JLD2.load(file)["confusions"]
conditions = JLD2.load(file)["conditions"]
operators  = JLD2.load(file)["operators"] 


##
function Base.unique(ctg::CategoricalArray)
    l = levels(ctg)
    newctg = CategoricalArray(l)
    levels!(newctg, l)
end
operators_labels = collect(map(x->String(x), keys(operators)))
models_labels = string.(1:4)
ctg = CategoricalArray(repeat(models_labels, inner = 7))
nam = CategoricalArray(repeat(operators_labels, outer = 4))
levels!(ctg, models_labels)
levels!(nam, operators_labels)
##

function get_chance_level()
    trials = round(Int, 50000)
    scores = falses(length(operators), trials)
    for x in 1:trials
        for o in 1:7
            res = rand([false,true])
            exp = rand(conditions) ∈ getfield(operators,keys(operators)[o])
            scores[ o,x] = exp
        end
    end
    return scores
end


function get_prediction_matrix(confusion_values)
    predictions = zeros(4,7)
    mapping = Dict(conditions[n]=>n for n in  1:4)
    for o in 1:7
        inputs, preds, target = confusion_values[o]
        counter = zeros(Int,4)
        for (i, p, t) in zip(inputs, preds, target)
            predictions[mapping[i],o] += (p +1)÷2
            counter[mapping[i]] +=1.
        end
        predictions[:,o] ./= counter
    end
    return predictions
end


scores=zeros(7,4)
bias=zeros(7,4,4)
predictions=zeros(7,4,4)
_scores = [scores_mem[model] for model in TN.labels]
_confusion = [confusions[model] for model in TN.labels]
expected = zeros(7,4)

for (n, op) in enumerate(operators)
    @show op
    for (m,cond) in enumerate(conditions)
        if cond in op
            expected[n,m] = 1
        end
    end
end
for m in 1:4
    for op in 1:4
        scores[:,m] .= _scores[m]
        p = get_prediction_matrix(_confusion[m])
        predictions[:,m,:] .= p'

    end
end

##


colors=[BLU BLU RED :black]
chance = reshape(repeat(mean(get_chance_level(), dims=2)[:,1],4),7,4)

scores
k_scores = zeros(size(scores))
for i in eachindex(scores)
    k_scores[i] = (scores[i] - chance[i])/(1-chance[i])
end

# p1= groupedbar(nam, chance .-0.5, group = ctg, xlabel = "Operators",        title = "", ylabel="Score", bar_width = 0.67, lw =0.5,background_legend=:white,legend=false, frame = :axes, ylims=(0.0,0.5), legendtitle="Tripod configuration", legendfontsize=12,  lc=:match, c=:grey,alpha=0.7, la=0.2, xticks=(1:7, operators_labels))

p1= groupedbar(nam, k_scores[:,:], group = ctg, xlabel = "Operators", ylabel="Score", bar_width = 0.67, lw =0.5,background_legend=:white,legend=false , frame = :axes, ylims=(0.0,1), legendtitle="Tripod configuration",       lc=:black, c=colors, xticks=(0.5:6.5, operators_labels));

plot!(p1, xlabel="Operators", ylabel = "Readout "*L"\kappa"*"-score", size=(600,260));
savefig(p1, plotsdir("dendritic_memory", "Fig7A.pdf"))


##
abIMMI = [1,2,3,4,5,6,7]
plots= [heatmap(predictions[:,x,:]', c=:amp) for x in [1,4,3]]
exp_mat = heatmap( expected', c=:amp, colorbar=false);
push!(plots, exp_mat)
p2 = plot(plots..., clims=(0,1), cbar=false, yticks=(1:4, string.(conditions)), tickfontsize=15, xticks=(1:7, string.(operators_labels[abIMMI])), xrotation=-90 );
savefig(p2, plotsdir("dendritic_memory", "Fig7B.pdf"))
p1
p2