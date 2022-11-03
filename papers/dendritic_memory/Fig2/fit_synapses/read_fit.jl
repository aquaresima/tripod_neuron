using Bootstrap
using Serialization

fit_soma() = x->()
fit_dend() = x->()
soma = deserialize(joinpath(@__DIR__, "fit", "soma.so"))
dendrite = deserialize(joinpath(@__DIR__, "fit", "dendrite.so"))


function print_synapses(p, label)
    println(label, " parameters are:")
    println("τrise: ", p[1] )
    println("τdecay: ", p[2] )
    println("gsyn: ", p[3]/TN.norm_synapse(p[1],p[2]) )
end

print_synapses(dendrite.t0, "dendrite")
print_synapses(soma.t0, "soma")
print_synapses([0.100, 15.0465,0.5652], "soma")
# τr = 5.468232581517444
# τd = 25.945900413159777
# gsyn = 0.21022818698969572
println(soma.t0)
#(15.046505347839819, 0.10000014137387765, 0.5652862549456903)

soma
