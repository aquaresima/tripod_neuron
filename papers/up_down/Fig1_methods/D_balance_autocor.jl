using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Statistics, StatsBase, StatsPlots, JLD2, Plots, RollingFunctions
include(projectdir("scripts", "plots", "default_plots.jl"))
using MultivariateStats
using Plots



simtime = 50_000
##
rate_in = 20
model = TN.models[3]
cond = TN.get_balance_conditions(model.ds..., rate_in; nmda = true, istdp_syn = true)
inputs = TN.make_spikes(simtime, 0; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)

r = TN.get_spike_rate(voltage);
c = TN.get_cv(voltage[1, 10000:end]);
@info "Model: $(model.ds) Hz, CV: $c, rate: $r";

plot(voltage[1, end-15_000:end], label = "")
plot!(voltage[2, end-15_000:end])
plot!(voltage[3, end-15_000:end])
plot!(voltage[1, end-15_000:end], c = :black, legend = false)


##

## Run a single example for distal distal model

model = TN.models[1]
cond = TN.get_balance_conditions(model.ds..., 30; nmda = true, istdp_syn = false)
inputs = TN.make_spikes(simtime; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
print(sum(inputs, dims = 2))
print(cond.synapses)
plot(voltage[1, :], title = mean(voltage[1, 1000:end]))
r = TN.get_spike_rate(voltage)
c = TN.get_cv(voltage)



##
simtime = 50_000
using Plots
model = TN.models[2]
β = 0
cond = TN.get_balance_conditions(model.ds..., 20; nmda = true)
inputs = TN.make_spikes(simtime, β; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
voltage
plot(
    autocor(TN.cap_voltage(voltage[1, 50000:end]), 0:10:1000),
    c = mpi_palette[1],
    label = "proximal NMDA",
)
cond = TN.get_balance_conditions(model.ds..., 20; nmda = false)
inputs = TN.make_spikes(simtime, 200; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
plot!(
    autocor(TN.cap_voltage(voltage[1, 50000:end]), 0:10:1000),
    c = mpi_palette[1],
    ls = :dash,
    label = "proximal no NMDA",
)


model = TN.models[1]
cond = TN.get_balance_conditions(model.ds..., 20; nmda = true)
inputs = TN.make_spikes(simtime, 200; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
voltage
plot!(
    autocor(TN.cap_voltage(voltage[1, 50000:end]), 0:10:1000),
    c = mpi_palette[2],
    label = "distal NMDA",
)
# plot!(CCA(cap_voltage(voltage[1,:]), cap_voltage([1,:]), 0:10:10000), c=mpi_palette[2], label="distal NMDA")
cond = TN.get_balance_conditions(model.ds..., 20; nmda = false)
inputs = TN.make_spikes(simtime, β; cond...)
voltage = TN.run_tripod(inputs, simtime; model..., do_spikes = true, cond...)
plot!(
    autocor(TN.cap_voltage(voltage[1, 50000:end]), 0:10:1000),
    c = mpi_palette[2],
    ls = :dash,
    label = "distal no NMDA",
)
# plot!(yscale=:log10)
plot!(legend = :topright, xlabel = "Time (ms)", ylabel = "Autocorrelation")# title="Autocorrelation of voltage")
plot!(ylims = (0.01, 1))#, yscale=:log10)
##
plot!(xlims = (0, 400), legend = :topright)

# myCCA=MultivariateStats.CCA(rand(1000),rand(1000),rand(1000,10), rand(1000,10), rand(10))
# MultivariateStats.fit(CCA,rand(1000), rand(1000))0

## Test input correlations

simtime = 1000_000
dt = TN.dt
τ = 1 / exp(-dt / 0.1)
noise = 0.0
inputs = zeros(round(Int, simtime))
for t = 1:round(Int, simtime)
    re = rand() - 0.5
    noise = (noise - re) / τ + re
    inputs[t] = maximum([0, noise])
    # inputs[t] =noise
end
inputs
# plot(inputs)
plot!(autocor(inputs, 1:10:1000))


std(inputs)

autocor(inputs, [1])
1 / τ
mean(inputs)
