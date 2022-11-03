using DrWatson
@quickactivate "Tripod"
using TripodNeuron
include(projectdir("scripts","dendritic_memory","default_plots.jl"))
include("epsp_dend_interaction.jl")
using StatsPlots

models = TN.models[[1,2]]

@. linear_model(x,p) = p[1]*x + p[2]

function gg_interaction(a,b; kwargs...)
	gg = gg_simulate(a, b;kwargs...)
	g1 = gg_simulate(a,0.;kwargs...)
	g2 = gg_simulate(0.,b;kwargs...)
	return (gg-(abs(g1)+abs(g2)))
end

function gege_interaction(model)
	function add_plot(func::Function, same::Bool, my_plot=nothing; label="", NMDA, model, max_y=1., kwargs...)
		gs = rand(1.:35,800,2)
		Δ, g = gg_synaptic_efficacy(gs=gs, same=same, func=func, NMDA=NMDA, model=model )
		fit = curve_fit(linear_model, g,Δ, [1.,0.])
		# plot!(g, g .* fit.param[1] .+ fit.param[2], ls=:dash, ; kwargs...)
		if my_plot == nothing
		else
			# scatter!(my_plot, g, Δ, mscolor=:black; kwargs...)
		end
		n = length(my_plot.series_list)

		return my_plot, g, fit 
	end
	p = plot()
	## Same branches
	p, g, fit2 = add_plot(gg_interaction, true, p; model=model, markershape=:cross, msc=BLU, markersize=8,markerstrokecolor=BLU, linewidth=2, color= BLU, label=L"g_eg_e", NMDA=false)

	p, g,  fit1 = add_plot(gg_interaction, true, p; model=model, markershape=:cross, msc=RED, markersize=8,markerstrokecolor=RED, linewidth=2, color= RED, label=L"g_eg_e (NMDA)", NMDA=true)

	## Different branches
	p, g, fit4 = add_plot(gg_interaction, false,p; model=model, markershape=:circle, msc=BLU, markersize=8,markerstrokecolor=BLU, linewidth=2, color= BLU, label=L"g_eg_e", NMDA=false)

	p, g, fit3 = add_plot(gg_interaction, false, p;model=model, markershape=:circle, msc=RED, markersize=8,markerstrokecolor=RED, linewidth=2, color= RED, label=L"g_eg_e (NMDA)", NMDA=true)

	return p,g, [fit1,fit2, fit3, fit4]
end

pdd, g, fitdd = gege_interaction(TN.H_distal_distal)
yticks!(-6:3:6, string.(-6:3:6))
ppp, g, fitpp = gege_interaction(TN.H_proximal_proximal)
yticks!(-6:3:6, string.(-6:3:6))

plot!(ppp, xlabel=L"|g_x g_x| (nS)^2", ylabel="", yaxis=true)
yticks!(-6:3:6, string.(-6:3:6))
plot!(pdd,ylabel=L"\Delta"*"EPSP (a.u.)")
plot!(pdd,tickfontsize=13, guidefontsize=18,background_color=:white)
p=plot(pdd,  ppp, layout=(1,2), legend=false,  markersize=5,tickfontsize=13, guidefontsize=18,background_color=:white, framestyle=:axes, xaxis=true)
plot(p, xticks=(0:200:800,string.(0:2:8)))
savefig(p,plotsdir("dendritic_memory","Fig5B_appendix.pdf"))


## make plots
function fit_goodness(g, fit)
	@assert(length(fit.resid) == length(g))
	mean(sqrt.((fit.resid./(g .* fit.param[1] .+ fit.param[2])).^2)/2)
end

function all_fit(g,myfit)
	σ =[]
	μ =[]
	for fit in myfit
		push!(σ, fit_goodness(g,fit))
		push!(μ, fit.param[1])
	end
	return σ, μ
end
σdd, μdd = all_fit(g,fitdd)
σpp, μpp = all_fit(g,fitpp)


# to get margin of error and confidence interval of each parameter at 5% significance level:


σs = reshape(hcat(σpp...,  σdd..., ),(4,2))
μs = reshape(hcat(μpp...,  μdd..., ),(4,2))
ctg = repeat(["AA \nNMDA", "AA \nAMPA", "AB \nNMDA", "AB \nAMPA"], outer = 2)
nam = repeat(["proximal" , "distal"] , inner = 4)


# Compute the reference residuals for a gaussian noise with Norm(x,x)

p1= groupedbar(ctg, 400μs, group = nam, xlabel = "Tripod configuration",
        title = "", ylabel=L"\Delta"*"EPSP (mV)", bar_width = 0.67,c = [ :brown :yellow], lw = 0, framestyle = :box, legend=false, ylims=(-2.,2.))

p2= groupedbar(ctg, σs, group = nam, xlabel = "Tripod configuration", ylabel = "Mean Square Residuals (σ)",
        title = "", bar_width = 0.67, c= [ :brown :yellow],
        lw = 0, framestyle = :box)
q = plot(p1,p2)
plot!(q,tickfontsize=13, guidefontsize=18,background_color=:white)

##
savefig(q,plotsdir("dendritic_memory","Fig5B.pdf"))

labels = ["AA NMDA", "AA AMPA", "AB NMDA", "AB AMPA"]
for (myfit,dend) in zip([fitdd, fitpp], ["Distal", "Proximal"])
	for (fit, label) in zip(myfit, labels)
	# We can estimate errors on the fit parameters,
	# to get standard error of each parameter:
		print("$dend dendrite: $label\n")
		# per_err = abs.(round.(sigma ./ fit.param *100, digits=2))[1]
		conf_level = 0.01
		value            =round.(fit.param[1], digits=5)
		sigma            =round.(stderror(fit)[1], digits=6)
		margin_of_error  =round.(margin_error(fit, conf_level)[1], digits=6)
		confidence_inter =round.(confidence_interval(fit, conf_level)[1], digits=6)
		residuals        =round(fit_goodness(g, fit), digits=3)
		print("Fit value $value \n")
		print("ΔEPSP with two 20nS co-active synapses:  $(400*value) \n mV")
		print("Standard error $sigma \n")
		print("margin of error: $margin_of_error\n")
		print("confidence interval: $confidence_inter\n")
		print("Residuals: $residuals")
		print("\n\n")
	end
end

savefig(q,plotsdir("dendritic_memory","Fig5B.pdf"))

