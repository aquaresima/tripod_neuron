# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"));

using Statistics

function input_protocol_tt(tt::Int64, comp::Int64 ;eff=1.f0, spikes::Matrix{Int})
    if spikes[comp,tt] >0
        return spikes[comp,tt]*eff
    else
        return 0
    end
end
using Random

Random.seed!(17)
simtime = 5_000
plots = []
for (label, model) in zip(TN.labels, TN.models)
    Random.seed!(17)
    spikes = TN.make_balance(simtime,model=model)
    spikes[2:3,:] .*=2.

    # spikes[1,:] = TN.PoissonInput(1000., simtime, TN.dt, neurons=1) 
    # voltage = TN.run_tripod(spikes, simtime; model..., AdEx=TN.AdExParams(u_r=-55))
    voltage = TN.run_tripod(spikes, simtime; model...)
    voltage = voltage[:,end-10_000:end]
    ν= round(TN.get_spike_rate(voltage), digits=1)
    p= plot(voltage[1,:], label="", title="$label \n rate: $ν Hz", lw=3)
    plot!(xticks=:none, ylims=(-70,25), top_margin=5Plots.mm)
    length(plots) ==0 && annotate!(-1000, 90, Plots.text("B", 18, :black, :bold))
    length(plots) == 0 && plot!(label="", title=L"\textbf{Active}"*"\n$label \n rate: $ν Hz", lw=3)
    length(plots) >2 && (plot!(xlabel="Time (s)", xticks=(range(0, 10_000,5), range(0,1,5))))
    push!(plots,p)
    @show label ν
end

p = plot(plots... ,size=(600,600), margin=3mm, layout=(4,1))

# The balance condition is about setting the excitatory/inhibitory balance of the neuron
# We estimate the balance when the output firing rate is less than 1Hz, the input rate is always fixed to 5kHz for all the experiment that use the `make_balance`  
simtime = 5_000
plots = []
for (label, model) in zip(TN.labels, TN.models)
    spikes = TN.make_balance(simtime,model=model)
    voltage = TN.run_tripod(spikes, simtime; model...)
    voltage = voltage[:,end-10_000:end]
    length(plots) == 3 && plot!(ylabel="                        Membrane potential (mV)")
    ν= round(mean(voltage[1,:]), digits=1)
    q= plot(voltage[1,:], label="", title="$label \nAv. membrane: $ν mV", lw=3)
    length(plots) ==0 && annotate!(-1000, -36, Plots.text("A", 18, :black, :bold))
    length(plots) == 0 && plot!(label="", title=L"\textbf{Inactive}"*"\n$label \nAv. membrane: $ν mV", lw=3)
    plot!(xticks=:none, ylims=(-70,-50), top_margin=5Plots.mm)
    length(plots) >2 && (plot!(xlabel="Time (s)", xticks=(range(0, 10_000,5), range(0,1,5))))
    push!(plots,q)
    @show label ν
end

q = plot(plots...,size=(600,600), margin=3mm, layout=(4,1))
#

p = plot(q,p, size=(1000,800))
savefig(plotsdir("dendritic_memory","Revision_balanced.pdf"))

p


