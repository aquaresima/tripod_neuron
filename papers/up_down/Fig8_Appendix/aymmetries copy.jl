##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using StatsBase
using JLD2, ProgressBars
using EllipsisNotation

include(projectdir("scripts", "plots", "default_plots.jl"))

##
get_data(title, type = "data") =
    load(datadir("up_down", "asymmetries", "rate_β=0$title.jld2"))
AMPA = (; Dict(Symbol(k) => x for (k, x) in get_data("AMPA"))...)
NMDA = (; Dict(Symbol(k) => x for (k, x) in get_data("NMDA"))...)

@unpack νs, ds = TN
ampa_color = :darkred
nmda_color = :darkblue
##
v = NMDA.rate
plot(mean(v[1, 2:end, 2:end, :], dims = (1, 2))[1, 1, :])
plot!(mean(v[1, 22:end, 2:22, :], dims = (1, 2))[1, 1, :])
plot!(mean(v[1, 22:end, 22:end, :], dims = (1, 2))[1, 1, :])
plot!(mean(v[1, 2:22, 2:22, :], dims = (1, 2))[1, 1, :])

ds[21]
##
using Random
simtime = 50_000
dd = (150, 150)
istdp_syn = false
do_spikes = true
_rate = 20
nmda = false

function run_example(nmda, dd; _rate, β)
    # Random.seed!(17)
    syn_model = nmda ? TN.human_synapses : TN.ampa_equivalent
    model = (
        ds = dd,
        syn_model = syn_model,
        species = "H",
        soma_only = dd[1] == 0 ? true : false,
        do_spikes = do_spikes,
    )
    cond = TN.get_balance_conditions(dd..., _rate, istdp_syn = istdp_syn, nmda = nmda)
    inputs = TN.make_spikes(simtime, β; cond...)
    voltage = TN.run_tripod(dd, inputs, simtime; synapses = cond.synapses, model...)
    return voltage
end


TN.models[1]
plot()
##
simtime = 50_000
samples = 10
βs = [0, 200, 400, 600]
z = zeros(length(βs), 2, 3, 3, length(νs))
for b in eachindex(βs)
    β = βs[b]
    for (nmda, c, n) in zip([true, false], [nmda_color, ampa_color], [1, 2])
        for m in eachindex(TN.models[1:3])
            ds = TN.models[m].ds
            cors = zeros(3, length(νs))
            for _ = 1:samples
                Threads.@threads for ν in eachindex(νs)
                    # @info "Rate:", νs[ν]
                    voltage = run_example(nmda, ds; _rate = ν, β = β)
                    voltage = voltage[:, end-40_0000:end]
                    c1 = cor(voltage[1, :], voltage[2, :])
                    c2 = cor(voltage[2, :], voltage[3, :])
                    c3 = cor(voltage[1, :], voltage[3, :])
                    cors[:, ν] .+= [c1, c2, c3]
                    # @warn "Mean: ", ν, mean(voltage[1,:])
                end
            end
            cors ./= samples
            z[b, n, m, :, :] .= cors[:, :]
        end
    end
end
data = @strdict cors = z βs = βs νs = νs models = TN.models
save(datadir("up_down", "asymmetries", "isi_CV.jld2"), data)
