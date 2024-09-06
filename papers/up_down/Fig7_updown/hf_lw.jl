using Pkg
Pkg.activate("@v1.10")
using StatsPlots, UncertainData, RollingFunctions

using DrWatson
@quickactivate "Tripod"
using TripodNeuron

using Plots, Revise, Random
using LsqFit
using Statistics
using HypothesisTests
using EffectSizes
using FFTW

include(projectdir("scripts", "plots", "default_plots.jl"))
include(projectdir("scripts", "up_down", "updown_analysis.jl"))
##

path = mkpath(plotsdir("up-down", "Figures"))
@unpack νs = TN
AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
model = TN.models[3]
simtime = 1_000_0

for n in 0:3
    τ = Int(5^n)    
    βs = 0:100:1000
    rates = 5:35

    ###
    function compute_lf_hf(vv)
        vv[vv.>AdEx.θ] .= AdEx.θ
        signal = vv[1:50:end] .- mean(vv)
        fs = 200          # Sampling rate (Hz)
        t0 = 1/fs              # Start time 
        tmax = simtime/1000          # End time       
        t = t0:1/fs:tmax;  
        @assert length(t) == length(signal)
        F = fftshift(fft(signal))
        freqs = fftshift(fftfreq(length(t), fs))
        x_hf = findall(x-> x> 30 && x<90, freqs)
        x_lf = findall(x-> x> 0 && x< 10, freqs)
        f = abs.(F)
        return mean(f[x_lf]) / mean(f[x_hf])
    end

    function compute_lf(vv)
        vv[vv.>AdEx.θ] .= AdEx.θ
        signal = vv[1:50:end] .- mean(vv)
        fs = 200          # Sampling rate (Hz)
        t0 = 1/fs              # Start time 
        tmax = simtime/1000          # End time       
        t = t0:1/fs:tmax;  
        @assert length(t) == length(signal)
        F = fftshift(fft(signal))
        freqs = fftshift(fftfreq(length(t), fs))
        x_lf = findall(x-> x> 0 && x< 5, freqs)
        f = abs.(F)
        return mean(f[x_lf])/ mean(f) 
    end

    function get_lf(_rate, model, β, simtime; nmda, sample=10 )
        return map(1:sample) do x
            cond = TN.get_balance_conditions(model.ds..., _rate; nmda = nmda)
            inputs = TN.make_spikes(simtime, β; cond..., τ = τ)
            voltage =
                TN.run_tripod(inputs, simtime; model..., cond..., AdEx = AdEx, do_spikes = true)
            compute_lf_hf(voltage[1,:]), compute_lf(voltage[1,:])
        end
    end

    ###
    lfhf_nmda = zeros(length(βs), length(rates))
    lfhf_ampa = zeros(length(βs), length(rates))
    lf_nmda = zeros(length(βs), length(rates))
    lf_ampa = zeros(length(βs), length(rates))
    for x in eachindex(βs)
        Threads.@threads for y in eachindex(rates)
            β = βs[x]
            _rate = y
            @show β, νs[_rate]
            data_lf = get_lf(_rate, model, β, simtime; nmda=true, sample=10 )
            lfhf_nmda[x,y] = mean([x[1] for x in data_lf])
            lf_nmda[x,y] = mean(x[2] for x in data_lf)
            data_lf = get_lf(_rate, model, β, simtime; nmda=false, sample=10 )
            lfhf_ampa[x,y] = mean([x[1] for x in data_lf])
            lf_ampa[x,y] = mean(x[2] for x in data_lf)
        end
    end
    data = @strdict lfhf_nmda lfhf_ampa lf_nmda lf_ampa βs rates νs= νs[rates] τ
    save(datadir("up_down", "bistability", "lfhf_τ$τ.jld2"), data)
end


##
plot()
colors = palette(:tab10,length(0:3))
for (n, c) in zip(0:3,  colors)
    τ = Int(5^n)    
    data = load(datadir("up_down", "bistability", "lfhf_τ$τ.jld2"))
    @unpack lfhf_ampa, lfhf_nmda =  dict2ntuple(data)
    @unpack lf_ampa, lf_nmda =  dict2ntuple(data)
    # color = palette(:bluesreds,21)
    # for x in 1:21
    plot!(mean(lf_nmda, dims=2), ribbon=std(lfhf_ampa, dims=2), c=c, msc=c,label="", shape=:circle, ms=2)
    plot!(mean(lf_ampa, dims=2), ribbon=std(lfhf_ampa, dims=2), c=c, msc=c,label="τ $τ ms", shape=:circle, ms=2)
    # end
    # scatter!([-1],[-1], label="τ $τ", c=c, ms=2, msc=c)
    hline!([10], label="", ls=:dash)
    # vline!([10], label="", ls=:dash)
end
plot!(xscale=:log, legend=:outerright)
# plot!([0,35],[0,35],ls=:dash, c=:lightgrey, label="")
# plot!([0,35],[0,35],ls=:dash, c=:lightgrey, label="")
plot!(xlabel="LF AMPA",ylabel="LF NMDA", size=(500,400))
##