# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


using Plots, ColorSchemes, LaTeXStrings, Revise, Measures
default(palette = palette(:tab10), frame=:axes, guidefontsize=18, 
                    tickfontsize=13, grid=false, margins=5mm, legend_foreground_color=:transparent)

RED = RGBA(0.889, 0.436,0.278,1)
BLU = RGBA(0,0.605,0.978,1)
GREEN = RGBA{Float32}(0.2422242978521988,0.6432750931576304,0.30444865153411527,1.0)
PURPLE = RGBA{Float32}(0.7644401754934356,0.4441117794687767,0.8242975359232758,1.0)
BROWN =
RGBA{Float32}(0.675544,0.555662,0.0942343,1.0)
BLUGREEN = RGBA{Float32}(4.82118e-7,0.665759,0.680997,1.0)
GREY = RGBA{Float32}(.4,.4,.4,1.)
color_list = [BLU RED GREEN PURPLE BROWN BLUGREEN]
colormap(L,c1=BLU, c2=RED) = reshape( range(c1, stop=c2,length=L), 1, L );
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

Spiketimes = Vector{Vector{Float32}}
function raster_plot(pop::Spiketimes; ax=plot(), kwargs...)
	_x, _y = Float32[], Float32[]
    y0 = Int32[0]
    for n in eachindex(pop)
		for ft in pop[n]
			push!(_x,ft*1e-3)
			push!(_y,n)
		end
	end
    plt = scatter!(_x , _y, m = (2, :black), leg = :none,
                  xaxis=("Time (s)" ), yaxis = ("neuron",); kwargs...)
	return plt
end