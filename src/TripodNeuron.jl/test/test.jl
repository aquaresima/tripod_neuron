
function TripodTest(TN)
    TN.Tripod()

    # mini tests
    v = TN.run_simulation(TN.null_input)
    @assert(abs(mean(v[1][1, :]) - TN.AdEx.Er) < 1)

    v = TN.simulate_neuron(TN.SST(), TN.null_input)
    @assert(abs(mean(v[1][:]) - TN.LIF_sst.Er) < 1)

    v = TN.simulate_neuron(TN.PV(), TN.null_input)
    @assert(abs(mean(v[1][:]) - TN.LIF_pv.Er) < 1)

    print("Tripod model loaded successfully")

    function constant_firing()
        tripod = TN.Tripod()
        spiked = false
        currents = [0.0, 0.0, 0.0]
        v = Vector()
        for i = 1:1000
            TN.exc_spike!(tripod.s)
            spiked = TN.update_tripod!(tripod, currents, spiked)
            push!(v, tripod.s.v)
        end
        tripod = TN.SST()
        spiked = false
        v = Vector()
        for i = 1:1000
            TN.exc_spike!(tripod)
            spiked = TN.update_lif_sst!(tripod, spiked)
            push!(v, tripod.v)
        end
    end
    constant_firing()

end

TripodTest(TN)

function test_synapse(;
    step::Int64,
    dt::Float32,
    tripod = nothing,
    soma = nothing,
    somaspike = false,
)

    if step == 1
        if tripod != nothing
            if somaspike
                exc_spike(tripod.s, 1.0)
                inh_spike(tripod.s, 1.0)
            else
                exc_spike(tripod.d[1], 1.0)
                inh_spike(tripod.d[1], 1.0)
            end
        end
        if soma != nothing
            exc_spike(soma, 1.0)
            inh_spike(soma, 1.0)
        end
    end
end
