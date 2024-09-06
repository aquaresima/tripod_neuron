using HDF5
using Plots
using RollingFunctions
include("../../../TripodNeuron/base/colors.jl")

function load_data(file)
    h5open(joinpath(@__DIR__, file), "r") do fid
        return read(fid)
    end
end

data = "data_synaptic_integration55.h5"
balance = "inhibitory_balance55.h5"


kie = load_data(balance)
kie["kie"]["active"]
νs = mydata["freqs"]
labels = mydata["labels"]
active = mydata["active"]
passive = mydata["passive"]

# %% markdown
# In the following we show how the inhibition/excitation rate depends on the dendritic geometry
# %% markdown
# Now that the correct balance has been computed, as such the membrane is stable around -55 mV, we compute the variance of the membrane while increasing the input rate
# %% codecell
# # vs = average_free( νs, KIE, 100);
# gtot, σgtot, spikes, corr_time = average_nonfree( νs, KIE, 100);
# membrane = dendritic_membrane(νs, KIE);
# %% codecell
# %% markdown
# The KIE factor necessary to keep the membrane voltage stable on -50 mV is computed by iteration.
# It increases because with higher rates the NMDA channels activate and more inhibition is required.
# The reason it drops in the 'distal-distal' configuration are unclear
# %% codecell
KIE_active = active["kie"]
KIE_passive = passive["kie"]
p = plot(
    νs,
    KIE_active ./ νs,
    xlabel = "input rate (kHz)",
    ylabel = "Kie ratio",
    label = labels,
    xscale = :log10,
    linewidth = 3,
    title = "Excitatory-Inhibitory balance",
    color = color_list,
    legendtitle = "Configuration",
)
p = plot!(
    νs,
    KIE_passive ./ νs,
    ls = :dash,
    xlabel = "input rate (kHz)",
    ylabel = "Kie ratio",
    label = labels,
    xscale = :log10,
    linewidth = 3,
    color = color_list,
    title = "Excitatory-Inhibitory balance",
    labels = "",
    legendfontsize = 12,
    legend = :topleft,
    legendtitlefontsize = 13,
)
ylims!(-1, 1)
savefig(p, joinpath(@__DIR__, "figures", "EI_balance.pdf"))
p
# %% markdown
# Now compute the average firing rate of each model
#
# %% markdown
# The membrane statistics ensure the previous computation was correct: each simulation is run with the appropriate KIE ratio.
# The variance of the membrane potential is connected with the spike rate and the variance of the soma-only condition is consistent with the paper from Kuhn 2004.
# Hereby the free-membrane dynamics
# %% codecell

vs_active = active["voltage"]
vs_passive = passive["voltage"]
p1p = plot(
    νs,
    vs_passive[:, :, 1],
    ribbon = vs_passive[:, :, 2],
    linewidth = 3,
    xscale = :log10,
    fillalpha = 0.2,
    legendtitle = "Model",
    legend = :bottomright,
    label = labels,
    title = "Membrane statistics",
)
p1a = plot(
    νs,
    vs_active[:, :, 1],
    ribbon = vs_active[:, :, 2],
    linewidth = 3,
    xscale = :log10,
    fillalpha = 0.2,
    legendtitle = "Model",
    legend = :bottomright,
    label = labels,
    title = "Membrane statistics",
)
# xticks!(vs)
ylabel!("Membrane potential (mV)")
xlabel!("Input rate (kHz)")
p2p = plot(
    νs,
    vs_passive[:, :, 2],
    lw = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Membrane variance",
    legend = false,
)
p2a = plot(
    νs,
    vs_active[:, :, 2],
    lw = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Membrane variance",
    legend = false,
)
ylabel!("Standard deviation")
xlabel!("Input rate (kHz)")
pactive = plot(p1a, p2a)
ppassive = plot(p1p, p2p)
savefig("Membrane_stat.pdf")
pa
# %% markdown
KIE_passive
νs
# The spike rate increases with the input rate when there are proximal compartments. When there are two distal dendritic branches the presence of inhibition is enough to stabilize the dendritic compartments below the spike threshold.
# %% codecell

Ke_soma = active["ke_soma"]
@show labels
n = 19
m = 1
tm = models[m]()
voltage = TN.simulate_fast(
    tm,
    1000,
    stimulation,
    d1 = νs[n],
    d2 = νs[n],
    s = 0.0,
    Kie = KIE_active[n, m],
)
voltage = TN.simulate_fast(
    tm,
    10000,
    stimulation,
    d1 = νs[n],
    d2 = νs[n],
    s = 0.0,
    Kie = KIE_active[n, m],
)
voltage = TN.simulate_fast(tm, 1000, stimulation, s = νs[n], Kie = KIE_active[n, m])
voltage = TN.simulate_fast(tm, 10000, stimulation, s = νs[n], Kie = KIE_active[n, m])

voltage = TN.simulate_fast(tm, 10000, stimulation, s = 0.0, Kie = KIE_active[n, m])

mean(voltage[1, :])

plot(voltage[1, :])

vs = free_membrane(νs = νs, KIE = KIE_active, samples = 1, models = models);

plot(vs[.., 1], ribbon = vs[.., 2])


mydata
spikes_active = active["spikes"]
spikes_passive = passive["spikes"]
# spikes_active = rates
# labels = ["dd" "dp" "pp" "ss"]
ps1 = plot(
    νs,
    spikes_active[:, :, 1],
    linewidth = 3,
    legend = :topleft,
    legendtitle = "Model",
    xscale = :log10,
    label = labels,
    title = "Spike rate",
    ylabel = "Output rate (Hz)",
)
ps2 = plot(
    νs,
    spikes_active[:, :, 2],
    lw = 3,
    legend = false,
    xscale = :log10,
    label = labels,
    title = "ISI \nCoeff. of variation",
)
# xticks!(vs)
ylabel!("CV ")
xlabel!("Input rate (kHz)")
ps = plot(ps1, ps2)
savefig(ps, "spike_cv.pdf")
# %% markdown
# The membrane potential of soma and dendrites in the spiking condition
# %% codecell
membrane = active["voltage"]
m1 = plot(
    νs,
    membrane[:, :, 1],
    lw = 3,
    ribbon = membrane[:, :, 2],
    legend = :bottomright,
    legendtitle = "Model",
    xscale = :log10,
    label = labels,
    title = "Soma membrane",
)
# m2 = plot(νs, membrane[:,1:3,1,2], lw=3, ribbon= membrane[:,:,2],  legend=false, xscale=:log10,label=labels, title= "Dendrite A membrane")
# m3 = plot(νs, membrane[:,1:3,1,3], lw=3, ribbon= membrane[:,1:3,2,3], legend=false,  xscale=:log10,label=labels, title= "Dendrite B membrane")
plot(m2, ylabel = "Membrane potential (mV)")
pm = plot(m2, m3, m1)
savefig(pm, "Membrane.pdf")
# %% codecell
corr_time_active = active["corrs"]
corr_time_passive = passive["corrs"]
p1 = plot(
    νs,
    corr_time_active[:, :, 1],
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Soma",
    legend = false,
)
p2 = plot(
    νs,
    corr_time_active[:, 1:3, 2],
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Dendrite A",
    legend = false,
)
p3 = plot(
    νs,
    corr_time_active[:, 1:3, 3],
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Dendrite B",
    legend = false,
)
# xticks!(vs)
ylabel!(
    p1,
    "Membrane \n
AutoCorrelation Time (ms)",
)
xlabel!(p2, "Input rate (kHz)")
pc = plot(p1, p2, p3, layout = (1, 3))
savefig(pc, "autocorrelation.pdf")
# %% markdown
# The effective membrane integration time is given by the the somatic capacitance divided all the active conductance in the neuron. The dendrites act as a fixed potential if the somatic integration is much faster than the dendritic integration:
#
# $$\tau_{soma} << \tau_{dend}$$
#
# But dendritic timescale are much faster than the soma when the input are onto the branches. Therefore we consider the relative timescale:
#
# $$\tau_r = \frac{\tau_{soma} \cdot \tau_{dend}}{\tau_{soma} + \tau_{dend}}$$
#
# $\tau_r$ is the effective time by which a dendritic input is integrated in the neuron. The time necessary to affects the somatic potential depends on the somatic integration timescale $\tau_{mem}$
# %% codecell
labels = ["dd" "dp" "pp" "ss"]
gax = []
syn_tot = sum(gtot[:, 1:3, 2:3], dims = 3)[:, :, 1]
g_s = [x().d[1].pm.g_ax + x().d[2].pm.g_ax for x in [dd dp pp]]
g_l = [1 / x().d[1].pm.Rm + 1 / x().d[2].pm.Rm for x in [dd dp pp]]

C_s = [1 / x().d[1].pm.g_ax + 1 / x().d[2].pm.C⁻ for x in [dd dp pp]]
# C_s
τeff_d = C_s ./ (syn_tot .+ g_s .+ g_l)
τeff_s = TN.AdEx.C ./ (g_s .+ TN.AdEx.gl)

# %% codecell
τeff = (τeff_d .* τeff_s) ./ (τeff_d .+ τeff_s)
# %% codecell
#pyplot()
labels = ["dd" "dp" "pp" "ss"]
gax = []
using LaTeXStrings
p = plot(
    νs,
    TN.AdEx.C ./ gtot[:, :, 1],
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = L"\tau_{eff}" * " soma",
)
ylabel!("Integration timescale (ms)")
s = plot(
    νs,
    τeff,
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = L"\tau_{eff}" * " dendrite",
    legend = false,
)
# xticks!(vs)
xlabel!("input rate")
pτ = plot(p, s, layout = (1, 2))
savefig(pτ, "taueff.pdf")
# %% codecell

# %% markdown
# Now let's consider the excitatory/inhibitory balance at each denedrite.
# First we compute the average incoming currents during each the free membrane stimulation
#
# %% codecell

function currents_membrane(νs, KIE)
    current = zeros(size(νs)[1], 4, 5, 3) #var/mean model, rate
    σcurrent = zeros(size(νs)[1], 4, 5, 3) #var/mean model, rate
    somacurrent = zeros(size(νs)[1], 4, 4) #var/mean model, rate
    σsomacurrent = zeros(size(νs)[1], 4, 4) #var/mean model, rate

    indexs = length(νs)
    Threads.@threads for n = 1:indexs
        ν = νs[n]
        function G(v, i)
            if i == 2
                return TN.NMDA_nonlinear.(v)
            else
                return 1.0
            end
        end
        for (m, tm) in enumerate([dd, dp, pp, ss])
            tm = tm()
            if (m < 4)
                voltage, currents, synapses = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = ν,
                    d2 = ν,
                    s = 0.0,
                    Kie = KIE[n, m],
                )
                voltage, currents, synapses = TN.simulate_tripod(
                    tm,
                    SIMTIME,
                    stimulation,
                    record_synapses = true,
                    d1 = ν,
                    d2 = ν,
                    s = 0.0,
                    Kie = KIE[n, m],
                )
            else
                voltage, currents, synapses = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = 0.0,
                    d2 = 0.0,
                    s = 2ν,
                    Kie = KIE[n, m],
                )
                voltage, currents, synapses = TN.simulate_tripod(
                    tm,
                    SIMTIME,
                    stimulation,
                    record_synapses = true,
                    d1 = 0.0,
                    d2 = 0.0,
                    s = 2ν,
                    Kie = KIE[n, m],
                )
            end
            E_rev = [0.0 0.0 -75.0 -90]
            for i = 1:4
                current[n, m, i, :] .= mean(
                    synapses[i, :, :] .* transpose(G(voltage, 2)) .*
                    transpose(E_rev[i] .- voltage),
                    dims = 1,
                )[
                    1,
                    :,
                ]
            end
            current[n, m, 5, :] .=
                1 / TN.AdEx.Rm * mean(TN.AdEx.Er .- transpose(voltage), dims = 1)[1, :]
            σcurrent[n, m, 5, :] .=
                std(1 / TN.AdEx.Rm * (TN.AdEx.Er .- transpose(voltage)), dims = 1)[1, :]
            somacurrent[n, m, :] .= mean(currents[:, :], dims = 2)[:, 1]
            σsomacurrent[n, m, :] .= mean(currents[:, :], dims = 2)[:, 1]
        end
    end
    return current, σcurrent, somacurrent, σsomacurrent
end

# %% codecell
currents, σcurrents, somacurrent, σsomacurrent = currents_membrane(νs, KIE);
#current = νs, models, synapses, compartments

total_curr(i) = abs.(
    sum(currents[:, :, 1:2, i], dims = 3)[:, :, 1] ./ sum(currents[:, :, 3:4, i], dims = 3)
)[
    :,
    :,
    1,
]
sigma_curr(i) =
    abs.(sum(σcurrents[:, :, :, i] ./ currents[:, :, :, i], dims = 3))[:, :, 1] .*
    total_curr(i)
# %% markdown
# The total current on each compartments are given by the flux of incoming ions through the synaptic conductances. The ratio between I/E currents is in the plot for each compartment. In the soma-soma condition all the synapses are on the soma (purple) but there are no synapses in the other three conditions. The global current onto the soma compartment is considered in the next plot.
# %% codecell
labels = ["dd" "dp" "pp" "ss"]
p1 = plot(
    νs,
    total_curr(1),
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    ribbon = sigma_curr(1),
    label = labels,
    title = "Soma",
    legend = false,
)
ylabel!("Excitatory/Inhibitory \ncurrents ratio")
p2 = plot(
    νs,
    total_curr(2),
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    ribbon = sigma_curr(3),
    label = labels,
    title = "Dendrite A",
    legend = false,
)
xlabel!("Input rate (KHz)")
p3 = plot(
    νs,
    total_curr(3),
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    ribbon = sigma_curr(3),
    label = labels,
    title = "Dendrite B",
    legend = false,
)
# xticks!(vs)
xlabel!("input rate")
plot!(p1, legend = :topright, legendtitle = "model")
pei = plot(p1, p2, p3, layout = (1, 3))
savefig("excitatory-inhibitory.pdf")
# %% markdown
# Ratio between the depolirizing currents coming through the dendritic branches and the leak current of the soma membrane.
# %% codecell
r = plot(
    νs,
    abs.(somacurrent[:, 1:3, 1] ./ soma_leak_curr),
    linewidth = 3,
    fillalpha = 0.2,
    xscale = :log10,
    label = labels,
    title = "Soma current",
)
ylabel!("Depolarizing / Hyperpolarizing \n currents ratio")
xlabel!("Input rate (KHz)")
pr = plot!(r, legend = :topleft, legendtitle = "model")
savefig(pr, "depvship.pdf")
# %% codecell
tm = pp()
plots = []
ν = 6
titles = ["distal-distal" "distal-proximal" "proximal-proximal" "soma-only"]
for (n, tm) in enumerate([dd, dp, pp])
    tm = tm()
    voltage, currents, synapses = TN.simulate_tripod(
        tm,
        1000,
        stimulation,
        d1 = νs[ν],
        d2 = νs[ν],
        s = 0.0,
        Kie = KIE[ν, n],
    )
    voltage, currents, synapses = TN.simulate_tripod(
        tm,
        1000,
        stimulation,
        d1 = νs[ν],
        d2 = νs[ν],
        s = 0.0,
        Kie = KIE[ν, n],
    )
    push!(plots, plot(voltage[1, :], label = "", title = titles[n]))
end
tm = ss()
voltage, currents, synapses =
    TN.simulate_tripod(tm, 1000, stimulation, d1 = 0.0, d2 = 0.0, s = 0.0, Kie = KIE[ν, 4])
voltage, currents, synapses = TN.simulate_tripod(
    tm,
    1000,
    stimulation,
    d1 = 0.0,
    s = νs[ν],
    d2 = 0.0,
    Kie = KIE[ν, 4],
)
push!(plots, plot(voltage[1, :], label = "", title = titles[4]))
plot(plots..., xticks = false, yticks = false)
plot!(plots[2], ylabel = "Membrane    potential    (mV)")
plot!(plots[4], xlabel = "time (ms)")
psoma = plot(plots..., layout = (4, 1), ylims = (-75, -20), titlefontsize = 10)
savefig(psoma, "somadynamics.pdf")
# %% markdown
# Nothing, just old code
# %% codecell
function exc_post(; tripod::Union{Nothing,TN.Tripod} = nothing, step, dt, target_d)
    if step == 100
        if target_d > 1
            TN.exc_spike(tripod.d[target_d-1], 50.0)
        else
            TN.exc_spike(tripod.s, 50.0)
        end
    end
end

function inh_post(; tripod::Union{Nothing,TN.Tripod} = nothing, step, dt, target_d)
    if step == 100
        if target_d > 1
            TN.inh_spike(tripod.d[target_d-1], 50.0)
        else
            TN.inh_spike(tripod.s, 50.0)
        end
    end
end

function after_spike(νs, KIE, target_d)
    voltage = zeros(size(νs)[1], 4, 2, 500) #var/mean model, rate
    indexs = length(νs)
    Threads.@threads for n = 1:indexs
        ν = νs[n]
        for (m, Tm) in enumerate([dd, dp, pp, ss])
            if (m < 4)
                tm = Tm()
                _, _, _ = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = ν,
                    d2 = ν,
                    s = 0.0,
                    Kie = KIE[n, m],
                )
                tm.s.v = -50
                voltageexc, currents, synapses = TN.simulate_tripod(
                    tm,
                    50,
                    exc_post,
                    record_synapses = true,
                    target_d = target_d,
                )
                tm = Tm()
                _, _, _ = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = ν,
                    d2 = ν,
                    s = 0.0,
                    Kie = KIE[n, m],
                )
                tm.s.v = -50
                voltageinh, currents, synapses = TN.simulate_tripod(
                    tm,
                    50,
                    inh_post,
                    record_synapses = true,
                    target_d = target_d,
                )

            else
                tm = Tm()
                _, _, _ = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = 0.0,
                    d2 = 0.0,
                    s = 2ν,
                    Kie = KIE[n, m],
                )
                tm.s.v = -50
                voltageexc, currents, synapses = TN.simulate_tripod(
                    tm,
                    50,
                    exc_post,
                    record_synapses = true,
                    target_d = target_d,
                )
                tm = Tm()
                _, _, _ = TN.simulate_tripod(
                    tm,
                    1000,
                    stimulation,
                    d1 = 0.0,
                    d2 = 0.0,
                    s = 2ν,
                    Kie = KIE[n, m],
                )
                tm.s.v = -50
                voltageinh, currents, synapses = TN.simulate_tripod(
                    tm,
                    50,
                    inh_post,
                    record_synapses = true,
                    target_d = target_d,
                )

            end
            voltage[n, m, 1, :] .= voltageexc[target_d, :]
            voltage[n, m, 2, :] .= voltageinh[target_d, :]
        end
    end
    return voltage
end

function findlocalmaxima(signal)
    inds = Int[]
    if length(signal) > 1
        if signal[1] > signal[2]
            push!(inds, 1)
        end
        for i = 2:length(signal)-1
            if signal[i-1] < signal[i] > signal[i+1]
                push!(inds, i)
            end
        end
        if signal[end] > signal[end-1]
            push!(inds, length(signal))
        end
    end
    inds
end
voltage = after_spike(νs, KIE, 3);
voltage .= voltage ./ voltage[:, :, :, 100];
extremes = []
for s = 1:6
    print(findlocalmaxima(voltage[s, 1, 1, 100:end]))
end
