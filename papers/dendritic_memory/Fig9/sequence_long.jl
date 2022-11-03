using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

using Logging, ProgressBars, StatsBase, StatsPlots
using MLJLinearModels, MLJ
include("stimuli.jl")
println("Run models")

function _classify(; β, interval, τinh, samples=200, model, A=5, make_plot=false)

    X,Y = run_batch(; A=A, model=model, β=β, interval=interval,
                            τinh=τinh, batch_size=samples)

    train_std = StatsBase.fit(ZScoreTransform, X, dims = 1)
    StatsBase.transform!(train_std, X)
    # @show mean(X, dims=1) std(X, dims=1)

    logistic = LogisticRegression(0.5)
    y =4*(Y .-1)
    theta = MLJLinearModels.fit(logistic, X, y)
    colors= [BLU, RED];c = round.(Int,Y)
    scatter(X * theta[2:3] .+ theta[1], c=c)

    X2,Y = run_batch(; A=A, model=model, β=β, interval=interval,
                            τinh=τinh, batch_size=samples)

    StatsBase.transform!(train_std, X2)
    colors= [BLU, RED];c = round.(Int,Y)
    ##
    y2 = X2 * theta[2:3] .+ theta[1]
    y0 =2*(Y .-1).-1
    ##
    z =sum(sign.(y0) .== - sign.(y2))/length(y0)
    ## if theta is zero: failed training
    z = sum(theta) == 0 ? 0.5 : z
    ## Solve mislabeling by taking the best guess of the model
    accuracy = abs(z- 0.5) +0.5

    if !make_plot
        return accuracy
    else
        ## make a nice plot
        class1 = [x for x in eachindex(y2) if y0[x]<=0]
        class2 = [x for x in eachindex(y2) if y0[x]>0]
        p = violin([zeros(length(class1)), ones(length(class2))], [y2[class1], y2[class2]], labels=["Order 1" "Order 2"])
        p = dotplot!([zeros(length(class1)), ones(length(class2))], [y2[class1], y2[class2]], c=:black, labels="")
        return accuracy, p , theta
    end
end

a,p, theta= _classify(; β=5, interval=30, τinh=50, model=TN.H_distal_proximal, make_plot=true)
plot(p, title="$a")

##

function get_scores(;model, inputs)
    @unpack intervals, βs, τinh = inputs
    scores = Array{Float64,2}(undef,length(intervals), length(βs))
    @info "Run model: $model"
    for r in ProgressBar(eachindex(βs))
    Threads.@threads for i in eachindex(intervals)
            S = _classify(model=model, β = βs[r], interval=round(Int,intervals[i]), τinh = τinh)
            scores[i,r] = S
        end
    end
    return scores
end

inputs =(
    βs = 0:1:20,
    intervals = 10 .^range(0.5,3,length=20),
    τinh = 50
)

data = Dict()
for (model, label) in zip(TN.models, TN.labels)
    @time push!(data,label =>get_scores(;model=model, inputs=inputs))
end

file =datadir("dendritic_memory","sequence_data.jld2")
data = @strdict data inputs
safesave(file, data)


