
function get_balance_conditions(
    d1,
    d2,
    rate;
    nmda = true,
    istdp_syn = false,
    shift = 1,
    syn_model::Symbol = :none,
)
    _rates = νs
    _models = ds

    if d1 + d2 <= 2
        d1 = 1
        d2 = 1
        soma_only = true
    else
        soma_only = false
    end

    if syn_model == :none
        if d1 > length(_models)
            d1 = argmin(map(l -> abs(d1 - l), _models))
        end
        if d2 > length(_models)
            d2 = argmin(map(l -> abs(d2 - l), _models))
        end
        if rate > length(_rates)
            rate = argmin(map(r -> abs(rate - r * 1000), _rates))
        end
        if nmda
            kie = getfield(balance_kie_rate, :nmda)
            gsyn = getfield(balance_kie_gsyn, :nmda)
            syn_model = default_synapses
        else
            kie = getfield(balance_kie_rate, :ampa)
            gsyn = getfield(balance_kie_gsyn, :ampa)
            syn_model = ampa_equivalent
        end
    else
        if istdp_syn
            throw("istdp_syn is not implemented for soma synapses")
        end
        kie = getfield(balance_kie_soma, syn_model)
        syn_model = getfield(soma_syn_models, syn_model)
    end

    l1 = _models[d1]
    l2 = _models[d2]

    if !istdp_syn
        synapses = [1, 1, 1, 1, 1, 1]
        kie_rate = [0, kie[rate, d1], kie[rate, d2]] .* shift
    else
        synapses = [1, 1, 1, 1, gsyn[rate, d1], gsyn[rate, d2]]
        kie_rate = [0, 1, 1]
    end
    return (
        ds = (l1, l2),
        rate = _rates[rate] * 1000,
        ie_ratio = kie_rate,
        synapses = synapses,
        istdp_syn = istdp_syn,
        soma_only = soma_only,
        syn_model = syn_model,
    )
end

function get_fluctuation_size(γs; rate = 1000)
    c = zeros(length(γs))
    for (i, γ) in enumerate(γs)
        x = () -> make_rates(4000, γ, rate = rate) |> rs -> var(rs) / mean(rs) / mean(rs)
        # y = () -> TN.make_rates(1000,γ,rate=1000) |> rs -> mean(rs)
        c[i] = mean([x() for _ = 1:500])
        # _m[i] = mean([y() for _ in 1:200])
    end
    return c
end


function get_efficacies(rate_range, effective_range)
    efficacy = Array{Float32,2}(undef, length(rate_range), length(effective_range))
    for n in eachindex(effective_range)
        for i in eachindex(rate_range)
            efficacy[i, n] = effective_range[n] / rate_range[i]
        end
    end
    return efficacy
end

function set_synapses(neurons, receptor::String, value::Float32)
    if isa(neurons, Array{Tripod,1})
        for neuron in neurons
            for d in neuron.d
                set_gsyn(getfield(d.syn, Symbol(receptor)), value)
            end
        end
    elseif isa(neurons, Tripod)
        for d in neurons.d
            set_gsyn(getfield(d.syn, Symbol(receptor)), value)
        end
    end
end
