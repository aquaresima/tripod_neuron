include("../../../TripodNeuron/simplified_model/simulation.jl")
include("../../../plasticity/inhibitory_plasticity/tripod_vogels.jl")
include("../../../TripodNeuron.jl")

##


data = read(h5open(joinpath(@__DIR__, "data/inhibitory_balance_45.h5"), "r"), "data")
kies = read(h5open(joinpath(@__DIR__, "data/inhibitory_balance_45.h5"), "r"), "kies")
νs = read(h5open(joinpath(@__DIR__, "data/inhibitory_balance_45.h5"), "r"), "freqs")
models = read(h5open(joinpath(@__DIR__, "data/inhibitory_balance_45.h5"), "r"), "labels")
opt_kies = zeros(size(data)[1:2])
for n in eachindex(νs)
    for m in eachindex(models)
        opt_kies[n, m] = kies[argmin(rollmean(data[n, m, :], 4))]
    end
end

# heatmap(data[:,10,:], clims=(0,0.5))
##
plot(data[15, 10, :])
plot!(0.5runmean(data[15, 10, :], 4) .+ 0.5rollmean(data[15, 10, :], 4), xlims = (50, 100))
plot!(runmean(data[15, 10, :], 4), xlims = (50, 100))
##

# opt_kies[n,m] = kies[argmin(runmean(data[n,m,:],10))]
ws = load(joinpath(@__DIR__, "data/istdp_balance.jld"), "data")
##
spikes = zeros(length(νs), length(models))
sigma = zeros(length(νs), length(models))

for n in eachindex(νs)
    @show n
    ν = νs[n]
    Threads.@threads for m in eachindex(TN.models_length)
        d = TN.models_length[m]
        model = (d, d)
        inputs = [0.0, ν, ν, 0.0, opt_kies[n, m] * ν, opt_kies[n, m] * ν]
        voltage = run_tripod(model, inputs, 5000)
        spikes[n, m] = TN.get_spike_rate(voltage[1, 10000:end])
        voltage = run_tripod_nospike(model, inputs, 5000)
        sigma[n, m] = var(voltage[1, 10000:end])
    end
end

##
p = heatmap(
    νs,
    models,
    spikes',
    c = :amp,
    ylabel = "Dendritic length (μm)",
    xlabel = "Dendritic length (μm)",
    xscale = :log,
)
q = heatmap(νs, models, sqrt.(sigma'), c = :amp, xscale = :log)
l = @layout [
    grid(1, 2)
    a{0.4h}
]

p = plot(p, q, plot(frame = :none), layout = l)
savefig(p, joinpath(@__DIR__, "spike_sigma.pdf"))

##

data = zeros(length(TN.models_length), length(TN.models_length))

n = 30
ν = νs[n]
for a in eachindex(TN.models_length)
    @show a
    Threads.@threads for b in eachindex(TN.models_length)
        d1 = TN.models_length[a]
        d2 = TN.models_length[b]
        model = (d1, d2)
        # kie1,kie2 = findfirst(x->x==model[1],TN.models_length), findfirst(x->x==model[2], TN.models_length)
        inputs = [0.0, ν, ν, 0.0, opt_kies[n, a] * ν, opt_kies[n, b] * ν]
        # synapses = [0.,1., 1., 0., TN.models_kie_gsyn[n], TN.models_kie_gsyn[m]]
        voltage = run_tripod_nospike(model, inputs, 5000)#, synapses)
        data[a, b] = mean(voltage[1, 5000:end])
    end
end


# savefig(p, joinpath(@__DIR__,"potential_3d.pdf"))

##

data_vogels = zeros(length(TN.models_length), length(TN.models_length))

ws
n = 18
ν = νs[n]
for a in eachindex(TN.models_length)
    @show a
    Threads.@threads for b in eachindex(TN.models_length)
        d1 = TN.models_length[a]
        d2 = TN.models_length[b]
        model = (d1, d2)
        # kie1,kie2 = findfirst(x->x==model[1],TN.models_length), findfirst(x->x==model[2], TN.models_length)
        inputs = [0.0, ν, ν, 0.0, ν, ν]
        synapses = [0.0, 1, 1, 0.0, ws[3, a, n], ws[3, b, n]]
        voltage = run_tripod(model, inputs, 1000, synapses = synapses)
        data_vogels[a, b] = TN.get_spike_rate(voltage[1, 5000:end])
    end
end




##
#pyplot()
p = heatmap(
    models,
    models,
    data,
    label = false,
    camera = (-45, 30),
    grid = false,
    ylabel = "Dendritic length (μm)",
    xlabel = "Dendritic length (μm)",
    cbar = false,
)
q = heatmap(
    models,
    models,
    data_vogels,
    label = false,
    xlabel = "Dendritic length (μm)",
    yaxis = false,
    cbar = false,
)

l = @layout [
    grid(1, 2)
    a{0.4h}
]
# ,   plot(frame=:none),layout=l )
s = plot(p, q, plot(frame = :none), layout = l)
savefig(s, joinpath(@__DIR__, "spikes_3d.pdf"))
