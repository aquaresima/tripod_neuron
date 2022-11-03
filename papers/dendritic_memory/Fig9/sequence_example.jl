include("stimuli.jl")

# Plot rate / CV data for all the models
markers = [:rect, :star5, :diamond, :circle]
plot()
for (model, label, marker) in zip(TN.models, TN.labels, markers)
    X,Y= run_batch(; model=model, β=15, interval=250, τinh=50, batch_size=100)
    colors = map(x->x==1 ? RED : BLU, Y)
    scatter!(X[:,1], X[:,2], msc=colors, c=colors, msw=3, ms=6, shape=marker,
    label=label)
end
p = plot!(legend=:topleft, ylabel="Firing rate (Hz)", xlabel="ISI CV", size=(400,400))
p
savefig(p, plotsdir("dendritic_memory", "Fig9F.pdf"))
p
