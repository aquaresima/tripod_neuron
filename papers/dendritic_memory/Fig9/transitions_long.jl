"""
Compute the correlation between ensemble-averaged firing rate and variations in the localization of input. 
The `test_signal` function does it in steps:
    1. Generate a signal.
    2. Compute the firing rate of a population of independent neurons exposed to the same siganl but with different spike inputs (balanced excitatory inhibitory inputs)
    3. Compute the variation of the signal (switch events in the journal article wordings).
    4. Measure the correlation between the signal and the population firing rate.
    5. Repeat the process for different intervals and tripod models.
"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include("stimuli.jl")
include(projectdir("scripts","dendritic_memory","default_plots.jl"));
using Logging

function test_signal(simetime, models, intervals, sig_func::Function, shifts, max_shift)
    scores = zeros(length(models), length(intervals), length(shifts))
    for (m, model) in enumerate(models)
        @info "Run model: $model"
        Threads.@threads for n in 1:length(intervals)
            interval = intervals[n]
            signal_raw = sig_func(simtime, interval)
            output, _ = spikes_integrator(model,signal_raw, simtime, samples=200)
            signal = signal_integrator(signal_raw, interval = 1:1.:simtime)
            shift  = measure_correlation(signal, output, max_shift=max_shift);
            scores[m,n,:] = shift
        end
    end
    return scores
end

file =datadir("dendritic_memory","transitions_corr_data.jld2")

simtime = 20_000
models = TN.models
intervals = 10 .^range(1,3.,length=25) #Hertz
max_shift=200
shifts = -max_shift+1:5:max_shift

@info "Testing oscillating signal"
scores_oscillation = test_signal(simtime, models,intervals, oscillating_signal, shifts, max_shift)
@info "Testing exp_decay  signal"
scores_exp_decay = test_signal(simtime, models,intervals, exp_decay_signal, shifts, max_shift)

data = @strdict scores_oscillation scores_exp_decay intervals models shifts
safesave(file,data)




