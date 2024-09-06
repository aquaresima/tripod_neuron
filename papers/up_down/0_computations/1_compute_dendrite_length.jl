"""
Compute the spiking rate of the tripod in all dendritic configurations
Synaptic input only on dendrites, and innhibition with same input rate than the excitation
"""
##
using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using JLD2, Statistics
using Logging

@info "Julia threads: $(Threads.nthreads())"


k_EI = 1.25
function test_length(; d1::Int, d2::Int, ν::Float64)
    voltage = TN.run_tripod((d1, d2), [0.0, ν, ν, 0.0, k_EI^-1 * ν, k_EI^-1 * ν], 10_000)
    r = TN.get_spike_rate(voltage)
    c = TN.get_cv(voltage[1, :])
    v = mean(voltage, dims = 2)
    return [r c v...]
end


function map_dspace(; ds, rates)
    data = zeros(5, length(ds), length(ds), length(rates)) ## rate/cv/voltage/soma/d1/d2
    Threads.@threads for i in eachindex(rates)
        ν = rates[i]
        @show ν
        for (n, d1) in enumerate(ds)
            for (m, d2) in enumerate(ds)
                data[:, n, m, i] = test_length(d1 = d1, d2 = d2, ν = ν)
            end
        end
    end
    return data
end

out = test_length(d1 = 200, d2 = 200, ν = 50.0)
##
println("Dendrite length test")
println("#####################")


data = map_dspace(ds = TN.ds, rates = TN.νs)
ds = TN.ds
νs = TN.νs
data = @strdict data νs ds reset = TN.AdEx.u_r k_EI
file = datadir("up_down", "rates", "dendritic_length_ds.jld2")
mkpath(dirname(file))
safesave(file, data)

##
function map_kspace(; ds, ks, rates)
    data = zeros(5, length(ks), length(ds), length(rates)) ## rate/cv/voltage/soma/d1/d2
    Threads.@threads for i in eachindex(rates)
        ν = rates[i]
        @show ν
        for (n, k) in enumerate(ks)
            for (m, d) in enumerate(ds)
                data[:, n, m, i] = test_length(k = Float64(k), d = d, ν = ν)
            end
        end
    end
    return data
end

function test_length(; k::Float64, d::Int, ν::Float64)
    voltage = TN.run_tripod((d, d), [0.0, ν, ν, 0.0, k^-1 * ν, k^-1 * ν], 10_000)
    r = TN.get_spike_rate(voltage)
    c = TN.get_cv(voltage[1, :])
    v = mean(voltage, dims = 2)
    return [r c v...]
end

ks_EI = 0.5:0.25:5
data = map_kspace(ds = TN.ds, ks = ks_EI, rates = TN.νs)
ds = TN.ds
νs = TN.νs
data = @strdict data ks = ks_EI νs ds reset = TN.AdEx.u_r
file = datadir("up_down", "rates", "dendritic_length_ks.jld2")
mkpath(dirname(file))
safesave(file, data)
