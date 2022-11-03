include("stimuli.jl")
using Random 

# Set parametes
simtime = 1000
T = 100
β = 15
τinh = 50
A=2

seed =  13
model_plots = []
for (model,label, letter) in zip(TN.models, TN.labels, "ABCD")
    # Run simulations in order 1 and order 2 with same random seed
    order = 1
    Random.seed!(seed)
    spikes = TN.make_balance(simtime, model=model)
    signal = three_phase_signal(simtime, T)
    stimuli = get_stimuli(signal,states=3)
    apply_three_phase_stimuli!(stimuli, spikes, order)
    voltage, v_i = TN.run_tripod_sequence(spikes, simtime,A =A ,τinh=τinh, β=β ; model...)

    cv =  round(TN.get_cv(voltage), digits=2)
    rate = round(TN.get_spike_rate(voltage), digits=2)

    q = plot(voltage[1,:] , lw=3, c=:black, title="Model: $label\n ABϵ, ISICV: $cv, rate: $rate", label="Soma")
    annotate!((-200, 90, Plots.text(letter, 18, :black, :bold)))
    if !(sum(model.ds) ==0 )
        q = plot!(voltage[2,:] , lw=3, label="Dendrite 1")
        q = plot!(voltage[3,:] , lw=3, label="Dendrite 2")
    end
    plot!(xticks=:none, legendfontsize=10, ylims=(-80,70))

    order = 2
    Random.seed!(seed)
    spikes = TN.make_balance(simtime, model=model) 
    signal = three_phase_signal(simtime, T)
    stimuli = get_stimuli(signal,states=3)
    apply_three_phase_stimuli!(stimuli, spikes, order)
    voltage, v_i = TN.run_tripod_sequence(spikes, simtime, A = A,τinh=τinh, β=β ; model...)
    
    cv =  round(TN.get_cv(voltage), digits=2)
    rate = round(TN.get_spike_rate(voltage), digits=2)
    p = plot(voltage[1,:] , lw=3, c=:black, title=" BAϵ, ISICV: $cv, rate: $rate")
    p= plot!(ylabel="                            Membrane potential (mV)")
    if !(sum(model.ds) ==0 )
        p = plot!(voltage[2,:] , lw=3)
        p = plot!(voltage[3,:] , lw=3)
    end
    plot!(xticks=:none, legend=false)

    xticks=(range(0, simtime/TN.dt, 5), range(0,simtime, 5))
    plot!(xticks=xticks, xlabel="Time (ms)")
    layout=@layout[a{0.61h} 
                   b{0.38h}]
    plot(q,p, layout=layout, margin=5mm)
    p= plot!(size=(600,600))
    push!(model_plots,p)
end
p= plot(model_plots..., size=(1200,1200), yticks=(-60:20:20, -60:20:20))
savefig(p, plotsdir("dendritic_memory", "Fig9_appendix.pdf"))
#

p