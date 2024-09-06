using RollingFunctions
using StatsBase
include("bimodal_kernel.jl")


function APH(
    voltage,
    title,
    h = 3;
    ax = plot(),
    filter = true,
    v_range = collect(-70:-40),
    kwargs...,
)
    p = Plots.histogram!(
        ax,
        voltage[1, :],
        normalize = true,
        bins = range(-70, stop = -40, length = 100),
        color = :black,
        label = "",
        title = title,
        titlefontsize = 15,
    )
    s = maximum(map(x -> !isnan(x) ? x : 0, p.series_list[1][:y]))
    if filter
        kde = globalKDE(h, voltage[1, :], v_range = v_range)
        plot!(v_range, kde, label = "", lw = 2; kwargs..., c = :darkorange)
        rate = round(TN.get_spike_rate(voltage), digits = 1)
        # annotate!([(-35, s/1.5, text(string(rate)*"Hz", :bottom,:center,:black ))])
    end
    return p
end


using StatsPlots, ColorSchemes
function add_aph(ax, voltage, shift = 0; color, rate::String)
    h =
        StatsBase.fit(Histogram, voltage, -80:0.5:-40, closed = :left) |>
        x -> normalize(x, mode = :probability)
    p = groupedbar!(
        ax,
        h.edges[1][1:end-1],
        [h.weights shift * ones(size(h.weights))],
        bar_position = :stack,
        labels = [rate ""],
        c = [color :transparent],
        lc = :transparent,
    )

end
function rate_histograms(rates_, simtime, model, β, nmda, do_spikes = true)
    CList = palette(:berlin, length(rates_))
    ax = plot(
        background_legend = :transparent,
        fg_color_legend = :transparent,
        legend = :bottomleft,
        legendtitle = "Input rates",
        legendtitlefontsize = 13,
        yticks = false,
        yaxis = false,
        xaxis = "Membrane potential (mV)",
        guidefontsize = 15,
        tickfontsize = 11,
        grid = false,
    )
    AdEx = TN.AdExParams(idle = 0.1, up = 0.1)
    annotate!(
        ax,
        [(-38, 0.15 * length(rates_) / 2, text("Neuron firing rate", rotation = -90, 15))],
    )
    for n in eachindex(rates_)
        color = CList[n]
        rate = rates_[n]
        shift = 0.15 * (length(rates_) - n)
        model[3] = rate
        cond = TN.get_balance_conditions(model..., nmda = nmda)
        inputs = TN.make_spikes(
            simtime,
            β;
            rate = cond.rate,
            soma = cond.soma,
            dend = cond.dends,
            ie_ratio = cond.ie_ratio,
        )
        voltage = TN.run_tripod(
            cond.model,
            inputs,
            simtime,
            AdEx = AdEx,
            do_spikes = do_spikes,
            soma_only = cond.soma,
            syn_model = cond.syn_model,
        )
        # μ = round.(mean(voltage[1,5000:end], dims=2)[:,1],digits=1)
        output = round(TN.get_spike_rate(voltage[1, :]), digits = 1)
        add_aph(
            ax,
            voltage[1, :],
            shift,
            color = color,
            rate = "$(round(rate, digits=0)) kHz",
        )
        annotate!((-40, shift, text("$output Hz", :right, :bottom, 10)))
    end
    model_name = sum(model[1:2]) > 0 ? "$(TN.balance_models[model[1:2]])" : "soma"
    name =
        (simtime = simtime / 1000, model = model_name, β = β, nmda = nmda) |>
        name ->
            savename("histogram", name, "pdf") |>
            name -> joinpath(plotsdir("up_down", "histograms"), name)
    savefig(ax, name)
end


##
function my_diag(array, dim = 1)
    @info "Diagonalizing array"
    @info "Dimeansions (compartment, ds, ds, rates): $(size(array))"
    @assert size(array)[2] == size(array)[3]
    _νs = size(array)[end]
    _ds = size(array)[2]
    new_dims = collect(size(array))
    new_array = zeros(_ds, _νs)
    # @show size(new_arrVay)
    for n = 1:length(new_array[:, 1])
        new_array[n, :] = array[dim, n, n, :]
    end
    return new_array
end
