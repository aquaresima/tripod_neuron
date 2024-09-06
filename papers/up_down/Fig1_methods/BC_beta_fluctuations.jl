using DrWatson
using Plots, Revise, Random
@quickactivate "Tripod"
using TripodNeuron
using Statistics
include(projectdir("scripts", "plots", "default_plots.jl"))
include(projectdir("scripts", "up_down", "default_plots.jl"))

##
simtime = 1_0000
rate = 1.5
r = []
βs =  0:60:480
for β = βs
    push!(r, TN.make_rates(simtime,β,rate=rate, seed = 29))
end
##
all_spikes = []
for i = eachindex(βs)
    β = βs[i]
	r_exc = zeros(3,simtime)
	r_inh = zeros(3,simtime)
    for n in 1:3
        r_exc[n,:] = r[i]
        # TN.make_rates(simtime,β,rate=rate)
        r_inh[n,:] = TN.make_rates(simtime,0,rate=rate)
    end
	rates = r_exc, r_inh
    for _ = 1:35
        spikes = TN.make_spikes(simtime, β, rate = 1.5, rates=rates)[2, :]
        push!(all_spikes, (findall(spikes .> 0)) .* TN.dt * 1000)
        # cv = findall(spikes.>0)*TN.dt |> TN.get_isi |>TN.CV_isi
        # @info cv
        # @show spikes
    end
    for _ = 1:10
        spikes = TN.make_spikes(1_0000, β, rate = 0.5)[2, :]
        push!(all_spikes, [])
    end
end

all_spikes

Spiketimes = Vector{Vector{Float32}}
p = raster_plot(Spiketimes(all_spikes))
#

up_r=Vector{}()
for i = eachindex(βs)
    β = βs[i]
    # r[i] .= TN.make_rates(1_0000, β, rate = rate) ./ 100 .+ β / 2)
    push!(up_r, 10*r[i] .+ β * 3/4)  #./ 100 .+ β / 2))
end

up_r

colors = palette(range(mpi_palette[5], stop = mpi_palette[8], length = length(up_r)))
colorgrad = cgrad(range(mpi_palette[5], stop = mpi_palette[8], length = length(up_r)))
shifts = [10, 40, 70, 100, 130, 160, 190, 220, 250]
for i in eachindex(up_r)
    plot!(-9999:0, fill(shifts[i], 10_000), c = :black, ls = :dash, lw = 2.5)
    p = plot!(-9999:0, up_r[i], c = colors[i], yticks = :none, label = "")
end
vline!([0], c = :black, ls = :dot)
plot!(
    ylabel = "Input rate (kHz)",
    xticks = (-9999:5000:10_001, 0:5:20),
    xlims = (-9999, 10_000),
)
cc = colorbar(
    1,
    1.5,
    2,
    c = colorgrad,
    digits = 2,
    topmargin = 5Plots.mm,
    yrotation = 90,
    ylabel = "Fluctuation size (β)",
    guidefontsize = 12,
    tickfontsize = 15,
    labelfontsize = 18,
)
layout = @layout [b{0.05w} a{0.9w,1.0h}]
p = plot(cc, p, layout = layout)

path = mkpath(plotsdir("up_down", "up_down", "robustness"))
savefig(p, joinpath(path, "input_fluctuations.pdf"))
p
##
data = load(datadir("up_down", "inputs", "cv_inputs.bson"))[:data]

data.τs[34]
q1 = plot(
    1 .+ (data.βs ./ 1000),
    mean(data.cvs[34, :, 1:5], dims = 2),
    ribbon = std(data.cvs[34, :, 1:5], dims = 2),
    label = "",
    xlabel = "Fluctuation size (β)",
    ylabel = "ISI CV",
    rightmargin = 13Plots.mm,
    xlims = (1, 2),
    leftmargin = 5mm,
)
plot!([1, 2], [1, 2], ls = :dash, c = :black, lw = 2.5, label = "")
q2 = contour(
    βs[2:end],
    τs,
    cvs[:, 2:end, 1],
    title = "CV of input spike trains",
    xlabel = "Fluctuations (β)",
    ylabel = "τ (ms)",
    legend = :topleft,
    ylims = (5, 1000),
    yscale = :log10,
    xscale = :log10,
    fill = true,
    c = :amp,
    rightmargin = 10mm,
)
layout = @layout [a{0.4w} b{0.6w}]
q = plot(q1, q2, layout = layout)

##
layout = @layout[
    b{0.7h}
    a{0.3h}
]
z = plot(p, q, layout = layout, size = two_columns_size .* (0.7, 1))
savefig(z, plotsdir("up_down/Figures", "Fig1.pdf"))




# ##For poster
# r = []
# for β in 10:60:250
# 	push!(r,TN.make_rates(1_000, β, rate=1000, seed=29))
# end
# q3 = plot(r, c=collect(palette(:roma,5))', lw=3, labels=collect(0:60:240)', legendtitle="β",
# ylabel="Input rate (kHz)", xlabel="Time (s)", yticks=((500,1000,2000), (0.5,1,2)), ylims=(500,2500), xlims=(0,1500), xticks=(0:250:1000, 0:0.25:1), frame=:axes)

# # q3 = 
# # plot!(ylabel="Input β", xlabel="", xticks=:none )
# data = load(datadir("up_down","robustness","cv_inputs.bson"))[:data]

# q1 = plot(data.βs,mean(data.cvs[:,:]',dims=2), ribbon=std(data.cvs[:,:]', dims=2), label="", xlabel="Fluctuation size (β)", ylabel="ISI CV")
# layout = @layout [a{0.3h}
# 				  b{0.4h}	
# 				  b{0.3h}]
# q = plot(q1,q3, layout=layout)

# path = mkpath(plotsdir("up_down","up_down","robustness"))
# savefig(q, joinpath(path,"input_specifics_poster.pdf"))
##

residual =
    (1 .+ (data.βs ./ 1000) .- mean(data.cvs[5:7, :]', dims = 2)) ./
    (1 .+ (data.βs ./ 1000) .+ mean(data.cvs[1:5, :]', dims = 2))
scatter(residual .* 100)
