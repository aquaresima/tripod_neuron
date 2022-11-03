using Plots
using Images
using RollingFunctions
using HDF5


"""
Get EPSP from Miles 1996 paper
"""


## Function
function get_center_col(column)
    newcol = zeros(size(column))
    yy = findall(x->x, column)
    if isempty(yy)
        return 0
    end
    y = round(Int,mean(yy))
    newcol[y] = true
    return newcol
end

function get_center(column)
    newcol = zeros(size(column))
    yy = findall(x->x, column)
    if isempty(yy)
        return 0
    end
    return mean(yy)
end

function ynormalize(neurite, yscale, yunit)
    neurite =  - neurite ./yscale * yunit
    neurite .-= neurite[1]
    return neurite
end
## Soma
miles_data=Images.load(joinpath(@__DIR__,"paper_pics/minimal.png"))
plot(miles_data)
white = RGBA{Float64}(1.,1.,1.,1.0)
noises = [(125:155,65:100), (100:125,40:60), ( 25:50,70:100), (85:125,150:250)]
for noise in noises
    miles_data[noise[1], noise[2]].=white
end
# First, remove the title
xscale = [75:85,150:220]
yscale = [10:100,250:300]
scale_x = deepcopy(miles_data[xscale[1], xscale[2]])
scale_y = deepcopy(miles_data[yscale[1], yscale[2]])
miles_data[xscale[1], xscale[2]] .= white
miles_data[yscale[1], yscale[2]] .= white
plot(miles_data)
dendrite = deepcopy(miles_data[1:100,:])
soma = deepcopy(miles_data[100:end,:])
plot(scale_y)
plot(scale_x)
##

soma_mask = map( x->x!=white, soma)
y  = filter(x-> x!=0,map(x->get_center_col(x), eachcol(soma_mask)))
soma_final = hcat(y...)
heatmap(soma_mask)
soma_filt  = filter(x-> x!=0,map(x->get_center(x), eachcol(soma_mask)))[2:end]

dendrite_mask = map( x->x!=RGBA{Float64}(1.,1.,1.,1.0), dendrite)
dendrite_final  = hcat(filter(x-> x!=0,map(x->get_center(x), eachcol(dendrite_mask)))...)
dendrite_filt =  filter(x-> x!=0,map(x->get_center(x), eachcol(dendrite_mask)))
heatmap(dendrite_mask)

# scale
mask = map( x->x!=white, scale_x)
xscale = length(filter(x-> x!=0,map(x->get_center(x), eachcol(mask))))
mask = map( x->x!=white, scale_y)
yscale  = length(filter(x-> x!=0,map(x->get_center(x), eachrow(mask))))
heatmap(hcat(y...))

myunit = 20
ratio = myunit / xscale

## rescale
soma = ynormalize(soma_filt, yscale, 2)
dendrite = ynormalize(dendrite_filt, yscale,2)
shorter = minimum(length.([soma, dendrite]))
xabs = ratio:ratio:shorter*ratio
soma = soma[1:shorter]
dendrite = dendrite[1:shorter]
p = plot(xabs,runmean(soma, 5), ylabel="mV", xlabel="ms");
p = plot!(xabs,runmean(dendrite, 5), ylabel="mV", xlabel="ms");
savefig(p, joinpath(@__DIR__,"miles_minimal.pdf"))

file = joinpath(@__DIR__,"miles_data_minimal.h5")
isfile(file) && rm(file)
# soma[soma .>0] .= 0.
# dendrite[dendrite .>0] .= 0.
h5open(file,"w") do fid
    fid["x"] = collect(xabs)
    fid["soma"] = soma
    fid["dend"] = dendrite
end
