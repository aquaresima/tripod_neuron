"""
This file runs all the scripts in the `papers` folder and reproduces all the figures in the `plot` folder.



"""

using DrWatson
@quickactivate "Tripod"
using TripodNeuron
using Logging
include(projectdir("scripts","dendritic_memory","default_plots.jl"))



root = projectdir("papers","dendritic_memory")
plots = plotsdir("dendritic_memory")
isdir(plots) || mkpath(plots)
data = datadir("dendritic_memory")
isdir(data) || mkpath(data)
dirs = [d for d in readdir(root) if isdir(joinpath(root,d))]

for dir in dirs
    path = joinpath(root,dir)
    for f in readdir(path)
        if !endswith(f,"long.jl") && isfile(joinpath(path,f))
            @info "running:" f
            include(joinpath(path,f))
        end
    end
end
