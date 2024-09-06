#ritic compartment length, coupled with the optimal one. 
## I decided to not put it in the paper.

# ##
# function maxrates(data, ls, rs)
#     diagonal = zeros(size(data)[2:3]...)
#     for x in 1:size(data)[2]
#         diagonal[x,:] .= data[x,x,:]
#     end
#     _rs = length(rs)
#     max_rate = runmean.([map(x-> maximum(x), eachcol(data[:,:,r])) for r in 1:length(rs)],15)
#     diag_rate = runmean.([diagonal[:,r] for r in 1:length(rs)],5)
#     mean_rate = runmean.([map(x-> mean(x), eachcol(data[:,:,r])) for r in 1:length(rs)],15)
#     mean_var = runstd.([map(x-> mean(x), eachcol(data[:,:,r])) for r in 1:length(rs)],15)

#     pmax = plot(legend=false)
#     pdiag = plot(legend=false)

#     min_rs =20
#     CList = palette(range(mpi_palette[5], stop=mpi_palette[8],length=length(νs[min_rs:end])))
#     for (r,c) in zip(min_rs:length(νs),CList)
#         pmax   =  plot!(pmax,  ls, max_rate[r], c=c)
#         scatter!([ls[argmax(max_rate[r])]], [maximum(max_rate[r])], c=c, markerstrokecolor=c, markersize = 5)
#         pdiag  =  plot!(pdiag, ls, diag_rate[r], c=c)
#         scatter!([ls[argmax(diag_rate[r])]], [maximum(diag_rate[r])], c=c, markerstrokecolor=c, markersize = 5)
#     end
#     # plot!(pmean, legend=false, ylabel="Firing rate (Hz)", xlabel="Dendritic length "*L"(\mu"*"m)", frame=:axes)
#     # plot!(pmax, legend=false, ylabel="Firing rate (Hz)", xlabel="Dendritic length "*L"(\mu"*"m)", frame=:axes)
#     cmap = cgrad(range(mpi_palette[5], stop=mpi_palette[8],length=length(νs[min_rs:end])))
#     return pmax, cmap
#     # return pmax, pdiag, pmean, cmap
# end

# pmax, cmap= maxrates(rs[2:end,2:end,:], ds[2:end], νs)
# hm,xss = colorbar_data(values = round.(Int,TN.νs[10:10:end]))
# z = plot(pmax, ylims=(0,130), xlabel="Dendritic length "*L"(\mu"*"m)", ylabel="Spike rate (Hz)", legend=false, frame=:box)
# inset = (1, bbox(0.1, 0.75, 0.3, 0.05, :bottom, :right))
# heatmap!(hm', c=cmap, frame=:box, yticks=:none, cbar=false, inset=inset, subplot=2, xticks=xss,  title="Input rate (kHz)", guidefontsize=15, tickfontsize=15, margin=0Plots.mm)
# z = plot(plot(frame=:none),z)
# # two_columns_size .* (1,1.5)
# # layout = @layout [  a{0.3h}
# #                     b{0.7h}]
# # plot(z, pq, layout=layout, size=two_columns_size .* (1,2), guidefontsize=18)
# savefig(z, joinpath(@__DIR__,"maxrate.pdf"))



# ## These functions compute other properties, but are not relevant for the last version of the paper.
# function best_coupling(data, ls, rs)
#     _rs = length(rs)
#     max_rate = [map(x-> argmax(x), eachcol(data[1,:,:,r])) for r in 1:length(rs)]

#     CList = range(BLU, stop=RED,length=_rs)
#     p = plot()
#     for r in 1:_rs
#         # p  =  plot!(ls,runmean(ls[max_rate[r]],3), c=CList[r], markerstrokecolor=CList[r])
#         p  =  scatter!(ls,ls[max_rate[r]]./ls, c=CList[r], markerstrokecolor=CList[r])
#         p  =  plot!(ls,ls[max_rate[r]]./ls, c=CList[r], markerstrokecolor=CList[r])
#     end

#     plot!(p, legend=false, ylabel="Best dendrititc length ratio", xlabel="Dendritic length "*L"(\mu"*"m)")

#     color_bar = reshape(repeat(1:_rs,2),(_rs,2))'
#     colormap(L,c1=BLU, c2=RED) = range(c1, stop=c2,length=L)
#     heatmap!(
#         color_bar,
#         inset = (1, bbox(0.05, 0.85, 0.30, 0.05, :bottom, :right)),
#         ticks = nothing,
#         c = ColorGradient(colormap(_rs)),
#         subplot = 2,
#         title="Input rate(kHz)",
#         titlefontsize=13,
#         xticks=([1,_rs],round.(rs[[1,_rs]], digits=1)),
#         colorbar = false,
#         bg_inside = nothing,
#         xscale=:log,
#         xlims =(0.8,31.6),
#         frame = :grid

#     )
#     return p
# end

# function fixed_dendrite(data, ls, rs; fixed_length::Int)
#     models = data[:,fixed_length,:,:]
#     _rs = length(rs)
#     diag_rate = runmean.([models[1,:,r] for r in 1:length(rs)],5)
#     CList = range(BLU, stop=RED,length=_rs)
#     pdiag = plot()
#     for r in 1:_rs
#         pdiag  =  plot!(pdiag, ls, diag_rate[r], c=CList[r])
#     end
#     plot!(pdiag, legend=false, ylabel="Spike rate (Hz)", xlabel="Dendritic length "*L"(\mu"*"m)")
#     scatter!(pdiag, ls[argmax.(diag_rate)], maximum.(diag_rate), c=CList, markerstrokecolor=CList, markersize = 5)

#     color_bar = reshape(repeat(1:_rs,2),(_rs,2))'

#     colormap(L,c1=BLU, c2=RED) = range(c1, stop=c2,length=L)
#     heatmap!(
#         color_bar,
#         inset = (1, bbox(0.05, 0.85, 0.30, 0.05, :bottom, :right)),
#         ticks = nothing,
#         c = ColorGradient(colormap(_rs)),
#         subplot = 2,
#         title="Input rate(KHz)",
#         titlefontsize=13,
#         xticks=([1,_rs],round.(rs[[1,_rs]], digits=1)),
#         colorbar = false,
#         bg_inside = nothing
#     )
#     return pdiag
# end

# ##
# plot!(pmax, legendfontsize=12, bglegend=:transparent, fglegend=:transparent,grid=false, guidefontsize = 18, tickfontsize=13)
# plot!(pdiag, legendfontsize=12, bglegend=:transparent, fglegend=:transparent,grid=false, guidefontsize = 18, tickfontsize=13)
# # pbest = best_coupling(mydata, ls, rs)
# # plot!(xscale=:log)
# f(x) = 250/x
# plot!(100:10:500,f.(100:10:500), c=:black,ls = :dash)
# f(x) = 100/x
# plot!(100:10:500,f.(100:10:500), c=:black, ls=:dash)
# plot!(pbest, legendfontsize=12, bglegend=:transparent, fglegend=:transparent,grid=false, guidefontsize = 18, tickfontsize=13)

# Plots.savefig(p1,   joinpath(@__DIR__,"rate_heatmap.pdf"))
# Plots.savefig(p2,   joinpath(@__DIR__,"cv_heatmap.pdf"))
# Plots.savefig(pmean,joinpath(@__DIR__,"mean_lines.pdf"))
# Plots.savefig(pmax, joinpath(@__DIR__,"max_lines.pdf"))
# Plots.savefig(pdiag, joinpath(@__DIR__,"diag_lines.pdf"))
# Plots.savefig(pbest, joinpath(@__DIR__,"best_coupling.pdf"))


# proximal = 6
# long_proximal = 11
# short_distal = 21
# distal= 31
# pprox = fixed_dendrite(mydata, ls, rs, fixed_length=long_proximal)
# pprox = fixed_dendrite(mydata, ls, rs, fixed_length=proximal)
# pshort = fixed_dendrite(mydata, ls, rs, fixed_length=short_distal)
# pdistal = fixed_dendrite(mydata, ls, rs, fixed_length=distal)


# rs[20]

# p1 = heatmap(ls, ls, mydata[1,:,:,15], colorbar=true, c=:amp)
# plot!(p1,ls,150*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,ls,300*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,ls,400*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,150*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# plot!(p1,300*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# plot!(p1,400*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# plot!(legend=false, tickfontsize=13, ylabel="Dendritic length "*L"(\mu m)",  xlabel="Dendritic length "*L"(\mu m)", guidefontsize=18, clabel="Spike Rate (Hz)")
# savefig(p1,joinpath(@__DIR__,"1khzheatmap.pdf"))

# ##
# p1 = heatmap(ls, ls, mydata[2,:,:,15], colorbar=true, clims=(0,1), c=:amp);
# plot!(p1,ls,150*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,ls,300*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,ls,400*ones(length(ls)), ls=:dash, c=:black, lw=3);
# plot!(p1,150*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# plot!(p1,300*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# plot!(p1,400*ones(length(ls)),ls, ls=:dash, c=:black, lw=3);
# annotate!(p1,[(160,405,text("distal",:black,15, :bottom,:left ))])
# annotate!(p1,[(160,305,text("short distal",:black,15, :bottom,:left ))])
# annotate!(p1,[(160,155,text("proximal",:black,15, :bottom,:left ))])
# plot!(legend=false, tickfontsize=13, ylabel="Dendritic length "*L"(\mu m)",  xlabel="Dendritic length "*L"(\mu m)", guidefontsize=18, clabel="Spike Rate (Hz)")
# savefig(p1,joinpath(@__DIR__,"1khzheatmap_cv.pdf")
