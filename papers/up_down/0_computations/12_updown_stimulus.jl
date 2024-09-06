using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Revise, Random
using LsqFit
using Statistics
using ProgressBars
include(projectdir("scripts", "plots", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_analysis.jl"))

simtime = 50_0000+ 1000
τ = 50
@unpack νs = TN
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
_rate = 11;
@show νs[_rate];
##
βs = collect(0:50:500)
model = TN.models[3]
for n in [2]
    τ = 5^n
    data = Vector{Any}(undef, length(βs))
    @show τ
    Threads.@threads for b in eachindex(βs)
        β = βs[b]
        Random.seed!(133)
        samples = map(1:50) do x
            cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
            inputs = TN.make_spikes(simtime, β; cond..., τ = τ, seed=x)
            voltage =
                TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
            v = voltage[1, 1_0001:end]
            v
        end
        vv = vcat(samples...)
        NMDA = updown_analysis(vv)
        
        Random.seed!(133)
        samples = map(1:50) do x
            cond = TN.get_balance_conditions(model.ds..., _rate; nmda = false)
            inputs = TN.make_spikes(simtime, β; cond..., τ = τ, seed=x)
            voltage =
                TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
            v = voltage[1, 1_0001:end]
            v
        end
        vv = vcat(samples...)
        AMPA = updown_analysis(vv)
        data[b] = @strdict NMDA AMPA τ
    end
    data = @strdict condition = data βs
    save(datadir("up_down", "bistability", "stimulus_betas_$τ.jld2"), data)
end

##
β = 200
rates = 9:2:25
data = Vector{Any}(undef, length(rates))
Threads.@threads for _rate in eachindex(rates)
    model = TN.models[3]
    # Random.seed!(133)
    cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
    inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
    voltage_NMDA =
        TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
    # histogram(voltage[1,:], bins=-80:-50, )
    NMDA = updown_analysis(voltage_NMDA[1, :])

    Random.seed!(133)
    cond = TN.get_balance_conditions(model.ds..., _rate; nmda = false)
    inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
    voltage_AMPA =
        TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
    # histogram(voltage[1,:], bins=-80:-50, )
    AMPA = updown_analysis(voltage_AMPA[1, :])

    data[_rate] = @strdict NMDA AMPA
end
data = @strdict condition = data rates
save(datadir("up_down", "bistability", "stimulus_rate_$τ.jld2"), data)

