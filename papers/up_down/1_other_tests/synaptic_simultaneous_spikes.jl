using QuadGK
using Plots
include("../TripodNeuron.jl")

τs = 1 / TN.Esyn_dend.AMPA.τd⁻
E0 = 0
C = TN.AdEx.C
V0 = TN.AdEx.Er
τm = TN.AdEx.τm
##
Θ(x::Real) = x > 0.0 ? x : 0.0
gsyn(t::Real, T::Real) = exp(-t / τs) #+ Θ(t-T)*exp(-(t-T)/τs)

h(t::Real, T::Real) = (1 / τm + gsyn(t, T) / 50)
g(t::Real, T::Real) = V0 / τm
## Solving the first ODE
μ(t::Real, T::Real) = quadgk(x -> exp(h(x, T)), 0, t, rtol = 1e-8)[1]

μ(10, 11)

function V(t::Real, T::Real)
    integral, err = quadgk(x -> μ(x, T) * g(x, T), 0, t, rtol = 1e-8)
    return integral / μ(t, T)
end
##
EPSPs = zeros(200, 1)
for t = 1:1:200
    for T = 1:1:1
        EPSPs[t, T] = V(t / 10, T / 10)
    end
end
##

heatmap(EPSPs, ylabel = "time(ms)", xlabel = "delay(ms)")

argmax(EPSPs, dims = 2)
