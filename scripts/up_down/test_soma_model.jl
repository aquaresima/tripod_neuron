

AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
function run_example(nmda, dd; β, _rate)
    Random.seed!(21)
    syn_model = nmda ? TN.human_synapses : TN.ampa_equivalent
    model = (
        ds = dd,
        syn_model = syn_model,
        species = "H",
        soma_only = dd[1] == 0 ? true : false,
        do_spikes = true,
    )
    cond = TN.get_balance_conditions(dd..., _rate, istdp_syn = false, nmda = nmda)
    inputs = TN.make_spikes(simtime, β; cond...)
    voltage, syns = TN.run_tripod(
        dd,
        inputs,
        simtime;
        synapses = cond.synapses,
        model...,
        rec_syn = true,
        AdEx = AdEx,
    )
    return voltage, syns
end


##
simtime = 1_0000
nmda = true
dd = (0, 0)
v, s = run_example(nmda, dd; _rate = 20, β = 300)
plot(v[1, :])
@assert sum(s[2, :, 2]) > 0
nmda = false
dd = (0, 0)
v, s = run_example(nmda, dd; _rate = 20, β = 300)
plot!(v[1, :])
@assert sum(s[2, :, 2]) == 0

rs_true = []
for ν = 1:40
    nmda = true
    dd = (0, 0)
    v, s = run_example(nmda, dd; _rate = ν, β = 300)
    r = TN.get_spike_rate(v)
    push!(rs_true, r)
end

rs_false = []
for ν = 1:40
    nmda = false
    dd = (0, 0)
    v, s = run_example(nmda, dd; _rate = ν, β = 300)
    r = TN.get_spike_rate(v)
    push!(rs_false, r)
end


plot(rs_true, label = "NMDA")
plot!(rs_false, label = "AMPA")
##
simtime = 1_0000
nmda = false
dd = (100, 100)
v = run_example(nmda, dd; _rate = 40, β = 0)
plot(v[2, :])
plot!(v[1, :])
# @assert sum(s[2,:,2]) > 0 
