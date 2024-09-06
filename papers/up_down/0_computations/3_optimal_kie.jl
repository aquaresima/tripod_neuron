using DrWatson
@quickactivate "Tripod"
using JLD2, HDF5, Statistics
using Revise
using TripodNeuron
## Load data from previous computations

get_data(title, type = "data"; suffix = "") =
    load(datadir("up_down", "balance", "inhibitory_balance_$title$suffix.jld2"))[type]

name = "_highres"
for μmem in [-55.0]
    data_raw = (
        istdp_NMDA = get_data("istdp_NMDA", suffix = "_ur-55.0"),
        istdp_AMPA = get_data("istdp_AMPA", suffix = "_ur-55.0"),
        dend_NMDA = get_data("NMDA", suffix = "_μmem$(μmem)$name"),
        dend_AMPA = get_data("AMPA", suffix = "_μmem$(μmem)$name"),
        νs = get_data("AMPA", suffix = "_μmem$(μmem)$name", "νs"),
        kies = get_data("AMPA", "kies"),
        # soma_55 = get_data("55.0_soma")
    )

    @unpack kies, νs = data_raw
    @unpack ds = TN
    @info "Frequencies" νs
    ##

    # heatmap(kies, νs,  sqrt.(data.dend["data"][:,15,:]),clims=(0,2), yscale=:log)
    # Get optimal_kies for each length and soma
    function weighted_average(kies, weights, threshold = 2.0)
        x = 0.0
        w = 0.0
        for (xx, ww) in zip(kies, weights)
            if ww < threshold
                x += xx * 1 / ww
                w += 1 / ww
            end
        end
        return w == 0 ? 0 : x / w
    end

    function get_opt_kie(data)
        @unpack kies, νs = data_raw
        opt_kies = zeros(size(data)[1:end-1])
        @info size(νs)
        for n in eachindex(νs)
            if ndims(opt_kies) > 1
                for m in eachindex(ds)
                    opt_kies[n, m] = weighted_average(kies, data[n, m, :])
                end
            else
                opt_kies[n] = weighted_average(data.kies, data.minima[n, :])
            end
        end
        return opt_kies
    end


    @info size(data_raw.dend_NMDA)
    @info size(data_raw.dend_AMPA)
    # size(data_raw.istdp_AMPA)
    opt_kies = Dict(
        name => get_opt_kie(getfield(data_raw, name)) for name in (:dend_AMPA, :dend_NMDA)
    )
    push!(
        opt_kies,
        Dict(
            name => mean(getfield(data_raw, name)[2:3, :, :], dims = 1)[1, :, :] for
            name in [:istdp_NMDA, :istdp_AMPA]
        )...,
    )
    opt_kies = (; opt_kies...)
    min_NMDA = findfirst(x -> any(opt_kies.dend_NMDA[x, :] .> 0), 1:length(TN.νs))
    min_AMPA = findfirst(x -> any(opt_kies.dend_AMPA[x, :] .> 0), 1:length(TN.νs))
    data =
        @strdict opt_kies models = ds νs = νs min_AMPA = min_AMPA min_NMDA = min_NMDA name

    filename = "optimal_IE_μmem$(μmem)$name.jld2"
    file = datadir("up_down", "balance", filename)
    safesave(file, data)
    file = datadir(TN.balance_path, filename)
    safesave(file, data)
end
