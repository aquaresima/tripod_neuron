
conditions = [:a, :b, :ab, :none]
operators=(
    A =   [:a, :ab],
    B =   [:b, :ab],
    AND = [:ab],
    OR =  [:a, :b, :ab],
    XOR = [:a, :b],
    IM =  [:none, :a, :ab ],
    MI =  [:none, :b, :ab ],
)


function get_condition_spikes(simtime, condition, model, symbol_time)
    spikes = TN.make_balance(simtime, model= model)
    ST = round(Int,4*symbol_time/TN.dt)
    spikes[1,:] = TN.PoissonInput(1000., simtime, TN.dt, neurons=1) 
    if condition ==:ab
        spikes[2,end-ST:end] .*= 2.
        spikes[3,end-ST:end] .*= 2.
    elseif (condition ==:a)
        spikes[2,end-ST:end] .*= 2.
    elseif (condition ==:b)
        spikes[3,end-ST:end] .*= 2.
    end
    return spikes
end


function feats_from_data(data::Matrix{Float64}, symbol_time)
    nr_feat = 5
    ll = round(Int,symbol_time/TN.dt)
    h = ll ÷ nr_feat
    features = zeros(Float64, nr_feat*2, size(data)[2])
    for (n,sample) in enumerate(eachcol(data))
        # sample[1:500] .= runmean(sample[1:500],50)
        features[1:nr_feat,n]     = [mean(sample[1+h*x:(x+1)*h]) for x in 0:nr_feat-1]
        features[nr_feat+1:end,n] = [mean(sample[500+1+h*x:500+(x+1)*h]) for x in 0:nr_feat-1]

    end

    return features
end


"""
Generate a sequence of states from a simulation with random configurations
:a :b :none :ab
Each state is composed of SYMBOL_TIME*dt features, that is the membrane value
of the soma during the exposition to the stimulus.
"""

function generate_sequence(model, symbols, symbol_time, u_r=-70)
    if isa(symbols, Number)  ## large number of trials for testing
        n_symbols = symbols
        my_conditions = rand(conditions,symbols)
    else ## one trial for each condition
        n_symbols = length(symbols)
        my_conditions = symbols
    end
    membrane   = zeros(round(Int,symbol_time/TN.dt), round(Int, n_symbols))
    adaptation = zeros(round(Int,symbol_time/TN.dt), round(Int, n_symbols))
    inputs = Array{Symbol, 1}(undef,n_symbols)

    idle_intervals = 5
    simtime = symbol_time * idle_intervals
    active_interval = round(Int,1+ symbol_time*(idle_intervals-1)/TN.dt)
    adexp_no_reset = TN.AdExParams(u_r=u_r)
    for x in 1:n_symbols
        condition=my_conditions[x]
        my_inputs = get_condition_spikes(simtime, condition, model, symbol_time)
        voltage, w_adapt = TN.run_tripod(my_inputs,simtime; model...,adapt=true, AdEx=adexp_no_reset)
        membrane[:,x]   = voltage[1,active_interval:end]'
        adaptation[:,x] = w_adapt[ active_interval:end]'
        inputs[x] = condition
    end
    return vcat(membrane, adaptation), inputs
end
