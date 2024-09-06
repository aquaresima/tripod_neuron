"""
This file reproduces Fig 2A and 2B of the paper
It is necessary to run compute_dendritic_length.jl before it
"""

include(joinpath(@__DIR__, "AB_plot_dendritic_length.jl"));

zz1
## Read data
data = load(datadir("up_down", "rates", "dendritic_length_ks.jld2"))
rs = data["data"][1, ..]
cv = data["data"][2, ..]
membrane = data["data"][3:end, ..]
νs = data["νs"]
ds = data["ds"]
ks = data["ks"]
data
# @info "Loaded data inh_ratio: $(data["inh_ratio"]), Ereset: $(data["reset"])"


##
p = plot()
rs
CList = collect(palette(:blues, 19))
z = []
s = []
for r = 1:40
    s = []
    for i = 1:14

        d = argmax(rs[i, 2:end, r]) + 1
        if rs[i, d, r] > 0
            scatter!(
                [νs[r]],
                [ds[d]],
                label = "",
                color = CList[i],
                markersize = 3,
                markerstrokewidth = 0.5,
                markerstrokecolor = CList[i],
            )
            push!(s, ds[d])
        end
        # scatter!([ds[d]], [rs[i,d,r]], label="", color=CList[i], markersize=3, markerstrokewidth=0.5, markerstrokecolor=CList[i])
    end
    s = isempty(s) ? [0] : s
    push!(z, mean(s))
end
plot!(νs[13:end], z[13:end], ls = :dash, lw = 2, color = :black, label = "mean")

p = plot!(
    xlims = (0.4, 50),
    ylims = (50, 500),
    ylabel = "Optimal dend. length (μm)",
    xlabel = L"ν_{exc} (Hz)",
    size = (600, 400),
)
plot!(xticks = ([1.5, 5], [1.5, 5]), legend = :topleft)
plot!(xticks = ([νs[13], 1.5, 5, 15, 50], [round(νs[13], digits = 1), 1.5, 5, 15, 50]))


annotate!(p, 0.05, 540, text("D", :center, 20))
plot!(xscale = :log)
inset = (1, bbox(0.1, 0.75, 0.3, 0.081, :bottom, :right))
inset = (1, bbox(0.1, 0.45, 0.08, 0.3, :bottom, :right))
p = colorbar!(
    p,
    inset,
    c = :blues,
    ks[1],
    ks[end],
    2,
    digits = 2,
    title = L"k_{E/I}",
    titlefontsize = 15,
    margin = 0Plots.mm,
    horizontal = false,
)

savefig(p, plotsdir("up_down/Figures", "Fig2C.pdf"))
p

CList = collect(palette(:blues, 20))'
rs
q = plot(
    νs,
    mean(rs[:, 2:end, :], dims = 2)[:, 1, :]',
    c = CList,
    msc = CList,
    legend = false,
    xscale = :log,
    ylabel = "Firing rate (Hz)",
    xlabel = L"ν_{exc} (Hz)",
    size = (600, 400),
)
# q = plot(mean(rs[1:end, 2:end, :], dims=2)[:,1,:]', c=CList, msc=CList, legend=false,  ylabel="Firing rate (Hz)", xlabel=L"ν_{exc} (Hz)", size=(600,400))
plot!(xticks = ([0.1, 0.7, 1.5, 5, 15, 50], [0.1, 0.7, 1.5, 5, 15, 50]))
annotate!(q, 0.010, 68, text("C", :center, 20))
inset = (1, bbox(0.70, 0.45, 0.08, 0.3, :bottom, :right))
q = colorbar!(
    q,
    inset,
    c = :blues,
    ks[1],
    ks[end],
    2,
    digits = 1,
    title = L"k_{E/I}",
    titlefontsize = 15,
    margin = 0Plots.mm,
    xrotation = -45,
    horizontal = false,
)

zz2 = plot(
    q,
    p,
    size = two_columns_size .* (1, 0.7),
    layout = (1, 2),
    margin = 0mm,
    bottommargin = 5mm,
    leftmargin = 15mm,
    rightmargin = 15mm,
)

layout = @layout [
    a{0.65h}
    b{0.25h}
]
z = plot(zz1, zz2, layout = layout, size = two_columns_size .* (1, 1.8))
savefig(z, plotsdir("up_down/Figures", "Fig2.pdf"))
