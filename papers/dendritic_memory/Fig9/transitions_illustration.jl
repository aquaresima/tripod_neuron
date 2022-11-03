"""
Compute the correlation between ensemble-averaged firing rate and variations in the localization of input. 
The `test_signal` function does it in steps:
    1. Generate a signal.
    2. Compute the firing rate of a population of independent neurons exposed to the same siganl but with different spike inputs (balanced excitatory inhibitory inputs)
    3. Compute the variation of the signal (switch events in the journal article wordings).
    4. Measure the correlation between the signal and the population firing rate.
    5. Repeat the process for different intervals and tripod models.
"""


include("stimuli.jl")



simtime = 6000
T = 150
shifts = 0:10:50
x0 = 30_000 # make plots after 3 seconds 
qs,ps  = [], []

# Use the same signals for all models
signal_exp = exp_decay_signal(simtime,  T)
signal = oscillating_signal(simtime, T)

p0 = plot(signal[x0:end],title="Oscillatory signal", c=:black)
q0 = plot(signal_exp[x0:end], title="Exp decay signal", c=:black)
push!(qs, q0)
push!(ps, p0)

# spiketimes are in ms not in timesteps
x0 = x0÷10
Clist = collect(palette(:roma,4))
Clist[end] = RGBA(0,0,0,1)
for (model, label,c, ) in zip(TN.models, TN.labels, Clist)
    ## generate a balanced exc-inh input spike train
    spikes = TN.make_balance(simtime, model=model)
    ## modulate it with the signal
    stimuli = get_stimuli(signal_exp)
    apply_stimuli!(stimuli, spikes)
    ## Run the simulation
    voltage1 = TN.run_tripod(spikes,simtime; model...)
    ## Compute the average firing rate and return the spiketimes for the raster plot
    rate, spiketimes = spikes_integrator(model, signal_exp, simtime, samples=10)
    spiketimes = [x for x in spiketimes[1] if x>x0][1:5:end]
    # isempty(spiketimes) && (spiketimes=[0])
    q=  plot(rate[x0:end] , c=c)
    @show label, mean(rate)
    q = Plots.scatter!(spiketimes .-x0, ones(length(spiketimes[1])), c=:black, ms=4)
    push!(qs,q)

    ## oscillating_signal
    spikes = TN.make_balance(simtime, model=model)
    sum(spikes,dims=2)/1000
    stimuli = get_stimuli(signal)
    apply_stimuli!(stimuli, spikes)
    voltage2 = TN.run_tripod(spikes,simtime, do_spikes=true; model...)
    rate, spiketimes = spikes_integrator(model, signal, simtime, samples=10)
    spiketimes = [x for x in spiketimes[1] if x>x0][1:5:end]
    # isempty(spiketimes) && (spiketimes=[0])
    p=  plot(rate[x0:end] , c=c)
    p = Plots.scatter!(spiketimes .-x0, ones(length(spiketimes[1])), c=:black, ms=4)
    push!(ps,p)
end
ps[2]
q = plot(ps..., layout=(5,1), size=(300,800), frame=:none, legend=false)
p = plot(qs..., layout=(5,1), size=(300,800), frame=:none, legend=false)
z = plot(p,q, size=(600,800))
savefig(z, plotsdir("dendritic_memory","Fig9C_illustration.pdf"))
