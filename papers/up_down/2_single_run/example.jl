Θ(x::Float64) = x > 0.0 ? x : 0.0
function alpha_function(t; t0, τ)
    if abs(t - t0) / τ > 5
        return 0.0f0
    else
        return (t - t0) / τ * exp32(1 - (t - t0) / τ) * Θ(1.0 * (t - t0))
    end
end

function convolve(spiketime::Vector{Float32}; interval::AbstractRange, τ = 100)
    rate = zeros(Float32, length(interval))
    @inbounds for i in eachindex(interval)
        v = 0
        t = interval[i]
        @simd for t0 in spiketime
            @fastmath if abs(t - t0) / τ < 5
                v += alpha_function(t, t0 = t0, τ = τ)
            end
        end
        rate[i] = v
    end
    return rate
end
TN.get_spike_rate(voltage[1, :])

plotly()
plot(voltage[1, 58800:60000])

AdEx = TN.AdExParams(idle = 0.1, up = 0.1, u_r = -55)
function run_example(nmda, dd; β, _rate)
    # Random.seed!(17)
    syn_model = nmda ? TN.human_synapses : TN.ampa_equivalent
    model = (
        ds = dd,
        syn_model = syn_model,
        species = "H",
        soma_only = dd[1] == 0 ? true : false,
        do_spikes = do_spikes,
    )
    cond = TN.get_balance_conditions(dd..., _rate, istdp_syn = istdp_syn, nmda = nmda)
    inputs = TN.make_spikes(simtime, β; cond...)
    voltage =
        TN.run_tripod(dd, inputs, simtime; synapses = cond.synapses, model..., AdEx = AdEx)
    return voltage
end


plot()
##
dd = (150, 150)
βs = 0:20:300
rate = zeros(2, 16)
νs[15]
nmda = true
dd = (300, 150)
ν = 15
for b in eachindex(βs)
    β = βs[b]
    # @info "Rate:", νs[ν]
    voltage = run_example(nmda, dd; _rate = 15, β = β)

    rate[1, b] = TN.get_spike_rate(voltage)
    rate[2, b] = TN.get_cv(voltage)
    # @warn "Mean: ", ν, mean(voltage[1,:])
end

plot(βs, rate[2, :])


##
simtime = 100_000
β = 250
dd = (400, 400)
voltage = run_example(nmda, dd; _rate = 15, β = β)

histogram(voltage[1, :], bins = -70:-50, normed = true, label = "NMDA")
cor(voltage[3, :], voltage[2, :])
aph = APH(voltage, "", 2, yticks = :none)
# plot!(aph, ylims=(0,m+m/4))
##
plot(voltage[1, :])
spiketimes = TN.get_spike_times(voltage)

histogram(diff(spiketimes), bins = 0:10:1000, normed = true, label = "NMDA")

intervals = diff(spiketimes[10:end] ./ 10)
spiketimes

cv2 = []
lls = []
for x = 2:length(intervals)
    z = 2 * abs(intervals[x-1] - intervals[x]) / (intervals[x-1] + intervals[x])
    ll = (intervals[x-1] + intervals[x]) / 2
    push!(cv2, z)
    push!(lls, ll)
end

scatter(lls, cv2, xlims = (0, 100), label = "NMDA")
mean(cv2)

##
get_data(syn) = load(datadir("up_down", "bistability", "critical_window_$syn.jld2"))
NMDA = get_data("NMDA")["data"]
heatmap(NMDA["proximal-proximal"].critical[2, :, :])
plot(mean(NMDA["proximal-proximal"].critical[2, :, :], dims = 2))
c = argmax(mean(NMDA["proximal-proximal"].critical[2, :, :], dims = 2))
M = maximum(mean(NMDA["distal-distal"].critical[2, :, :], dims = 2))

NMDA["proximal-proximal"].critical[2, c[1], c[2]]
TN.νs[6]

hline!([6], label = "")
