"""
Compute the spiking rate of the tripod in all dendritic configurations
Synaptic input only on dendrites, and innhibition with same input rate than the excitation
"""
##
include("../../base.jl")
using HDF5


function test_ν(; d::Int, ν1::Float64, ν2::Float64)
    model = "H. 1->s, $d,4; 2->s, $d,4"
    # tripod = TN.Tripod(TN.TripodModel(model))
    voltage = TN.run_tripod((d, d), [0.0, ν1, ν1, 0.0, ν2, ν2], 10_000)
    r = TN.get_spike_rate(voltage)
    c = TN.get_cv(v = voltage[1, :])
    v = mean(voltage, dims = 2)
    return [r c v...]
end


function map_νspace(; ds, νs)
    data = zeros(5, length(νs), length(νs), length(ds)) ## rate/cv/voltage/soma/d1/d2
    Threads.@threads for i in eachindex(ds)
        d = ds[i]
        @show d
        for (n, ν1) in enumerate(νs)
            for (m, ν2) in enumerate(νs)
                data[:, n, m, i] = test_ν(ν1 = ν1, ν2 = ν2, d = d)
            end
        end
    end
    return data
end

out = test_ν(d = 200, ν1 = 10.0, ν2 = 7.0)
##
println("Dendrite length test")
println("#####################")
ds = 100:10:500
νs = exp.(range(log(0.1), log(50), length = 40))
data = map_νspace(ds = ds, νs = νs)

file = joinpath(@__DIR__, "data/dendritic_freq.h5")
isfile(file) && (rm(file))
h5open(file, "w") do fid
    fid["rate"] = data[1, :, :, :]
    fid["cv"] = data[2, :, :, :]
    fid["soma"] = data[3, :, :, :]
    fid["d1"] = data[4, :, :, :]
    fid["d2"] = data[5, :, :, :]
    fid["length_range"] = collect(ds)
    fid["rate_range"] = collect(νs)
end
