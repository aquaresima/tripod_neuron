include("../TripodNeuron.jl")
using Plots

a = TN.PoissonInput(1400, 1000, TN.dt)
b = TN.PoissonInput(1400, 1000, TN.dt)
c = TN.PoissonInput(1000, 1000, TN.dt)
d = TN.PoissonInput(1000, 1000, TN.dt)
function protocol(; tripod = nothing, step::Int64, dt::Float32)
    (a[1, step]) && (TN.exc_spike!(tripod.d[1]))
    (b[1, step]) && (TN.exc_spike!(tripod.d[2]))
    (c[1, step]) && (TN.inh_spike!(tripod.d[1]))
    (d[1, step]) && (TN.inh_spike!(tripod.d[2]))
end

function protocol_rnd(; tripod = nothing, step::Int64, dt::Float32)

    TN.set_rate!(tripod.d[1], 2.0)
    TN.set_rate!(tripod.d[1], -1.0)
    TN.set_rate!(tripod.d[2], 2.0)
    TN.set_rate!(tripod.d[2], -1.0)
end

##
rates = zeros(200, 3)
cvs = zeros(200, 3)

t = TN.Tripod(TN.H_distal_proximal)
for h = 1:200
    t = TN.Tripod(TN.H_distal_proximal)
    voltage = TN.simulate_fast(t, 1000, protocol_rnd)
    rates[h, 1] = TN.get_spike_rate(voltage)
    cvs[h, 1] = TN.get_cv(voltage)
end

t = TN.Tripod(TN.H_distal_proximal)
voltage = TN.simulate_fast(t, 1000, protocol_rnd)
plot(voltage[2, 1:1000])
plot!(voltage[3, 1:1000])

t = TN.Tripod(TN.H_distal_proximal)
@eval TN.AdEx.idle_dt = 20

for h = 1:200
    t = TN.Tripod(TN.H_distal_proximal)
    voltage = TN.simulate_fast(t, 1000, protocol_rnd)
    rates[h, 2] = TN.get_spike_rate(voltage)
    cvs[h, 2] = TN.get_cv(voltage)
end

t = TN.Tripod(TN.H_distal_proximal)
@eval TN.AdEx.up_dt = 10

for h = 1:200
    t = TN.Tripod(TN.H_distal_proximal)
    voltage = TN.simulate_fast(t, 1000, protocol_rnd)
    rates[h, 3] = TN.get_spike_rate(voltage)
    cvs[h, 3] = TN.get_cv(voltage)
end
cvs[isnan.(cvs)] .= 0.0

histogram(rates[:, 1], label = "spike only")
histogram!(rates[:, 2], alpha = 0.5, label = "ref period")
histogram!(rates[:, 3], alpha = 0.5, label = "bAP + ref period")
plot!(xlabel = "Spike rate (Hz)", ylabel = "Number of realizations")
savefig("spike_count.png")

##
histogram(cvs[:, 1], bins = 0:0.01:1, label = "spike only")
histogram!(cvs[:, 2], bins = 0:0.01:1, alpha = 0.5, label = "ref period")
histogram!(cvs[:, 3], bins = 0:0.01:1, alpha = 0.5, label = "bAP + ref period")

##
t = TN.Tripod(TN.H_distal_proximal)
@eval TN.AdEx.up_dt = 0
@eval TN.AdEx.idle_dt = 20
voltage = TN.simulate_fast(t, 1000, protocol)
plot(voltage[1, :], c = :black, lw = 3, label = "Soma")
plot!(voltage[2, :], c = TN.RED, lw = 2, label = "Distal")
p = plot!(voltage[3, :], c = TN.GREEN, lw = 2, label = "Proximal")
plot!(p, title = "ref period")

t = TN.Tripod(TN.H_distal_proximal)
@eval TN.AdEx.up_dt = 10
@eval TN.AdEx.idle_dt = 20
voltage = TN.simulate_fast(t, 1000, protocol)
plot(voltage[1, :], c = :black, lw = 3, label = "Soma")
plot!(voltage[2, :], c = TN.RED, lw = 2, label = "Distal")
q = plot!(voltage[3, :], c = TN.GREEN, lw = 2, label = "Proximal")
plot!(q, title = "bAP + ref period")

t = TN.Tripod(TN.H_distal_proximal)
@eval TN.AdEx.up_dt = 0
@eval TN.AdEx.idle_dt = 0
voltage = TN.simulate_fast(t, 1000, protocol)
plot(voltage[1, :], c = :black, lw = 3, label = "Soma")
plot!(voltage[2, :], c = TN.RED, lw = 2, label = "Distal")
s = plot!(voltage[3, :], c = TN.GREEN, lw = 2, label = "Proximal")
plot!(s, title = "spike only")

##

z = plot(q, p, s, frame = :box, layout = (3, 1), legend = false)
savefig(z, "membrane.png")


TN.get_spike_rate(voltage)
