using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsPlots, JLD2, Plots, RollingFunctions
include(projectdir("scripts", "up_down", "default_plots.jl"))
using Random
using Plots

##
# Random.seed!(1234)
simtime = 1_0000
model = TN.models[1]
_rate = 14
@info "Exc rate: $(TN.νs[15])"
cond = TN.get_balance_conditions(model.ds..., _rate; nmda = true)
inputs = TN.make_spikes(simtime, 0; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
plot(voltage[1, :], title = mean(voltage[1, 1000:end]))
# plot(voltage[1,end-80000:end], title= mean(voltage[1,1000:end]))
voltage2 = voltage
# TN.get_spike_rate(voltage)
#

plot()
# Random.seed!(1234)
cond = TN.get_balance_conditions(model.ds..., _rate; nmda = false, syn_model = :kuhn)
inputs = TN.make_spikes(simtime, 300; cond...)
voltage = TN.run_tripod(model.ds, inputs, simtime; do_spikes = true, cond...)
print(sum(inputs, dims = 2))
print(sum(inputs, dims = 2))
# print( cond.synapses)
# print( cond.synapses)
@show sum(voltage .- voltage2)
plot!(voltage[1, :], title = mean(voltage[1, 1000:end]))
# plot!(voltage[1,end-80000:end], title= mean(voltage[1,1000:end]))
TN.get_spike_rate(voltage)
# TN.get_cv(voltage)
# histogram(voltage[1,:], bins=100, title="Soma voltage")
##
# simtime = 500_000
# using Plots
# model = TN.models[2]
# β = 100
# cond=TN.get_balance_conditions(model.ds..., 20; nmda=true)
# inputs = TN.make_spikes(simtime, β;  cond...)
# voltage =  TN.run_tripod(inputs,simtime; model...,do_spikes=true, cond...)
# voltage
# plot(autocor(cap_voltage(voltage[1,:]), 0:10:10000), c=mpi_palette[1], label="proximal NMDA")
# cond=TN.get_balance_conditions(model.ds..., 20; nmda=false)
# inputs = TN.make_spikes(simtime, 200;  cond...)
# voltage =  TN.run_tripod(inputs,simtime; model...,do_spikes=true, cond...)
# plot!(autocor(cap_voltage(voltage[1,:]), 0:10:10000), c=mpi_palette[1], ls=:dash, label="proximal no NMDA")
# # plot!(yscale=:log10)

# model = TN.models[1]
# cond=TN.get_balance_conditions(model.ds..., 20; nmda=true)
# inputs = TN.make_spikes(simtime, 200;  cond...)
# voltage =  TN.run_tripod(inputs,simtime; model...,do_spikes=true, cond...)
# voltage
# plot!(autocor(cap_voltage(voltage[1,:]), 0:10:10000), c=mpi_palette[2], label="distal NMDA")
# cond=TN.get_balance_conditions(model.ds..., 20; nmda=false)
# inputs = TN.make_spikes(simtime, β;  cond...)
# voltage =  TN.run_tripod(inputs,simtime; model...,do_spikes=true, cond...)
# plot!(autocor(cap_voltage(voltage[1,:]), 0:10:10000), c=mpi_palette[2], ls=:dash, label="distal no NMDA")
# # plot!(yscale=:log10)
# plot!(legend=:topright, xlabel="Time (ms)", ylabel="Autocorrelation")# title="Autocorrelation of voltage")
# plot!(ylims=(0.01,1), yscale=:log10)
