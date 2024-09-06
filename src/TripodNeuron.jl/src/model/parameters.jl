# Copyright (c) 2022 Author Name (AQ)
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT


const dt = 0.1f0

include("units.jl")
include("synapses.jl")
include("dendrites.jl")
include("soma.jl")

## Set neuronal parameters
const AdExTripod = AdExParams(Er = -70.6, u_r = -70.6)#, up=1., idle=0)
const AdEx = AdExTripod
const postspike = PostSpike()

## Set physiology
"""
(human) Eyal, G.; Verhoog, M. B.; Testa-Silva, G.; Deitcher, Y.; Lodder, J. C.; Benavides-Piccione, R.; Morales, J.; DeFelipe, J.; de Kock, C. P.; Mansvelder, H. D.; Segev, I. Unique Membrane Properties and Enhanced Signal Processing in Human Neocortical Neurons. eLife 2016, 5, e16553. https://doi.org/10.7554/eLife.16553.
"""
const HUMAN = Physiology(200Ω * cm, 38907Ω * cm^2, 0.5μF / cm^2)

"""
(mouse) (1) Dasika, V. K.; White, J. A.; Colburn, H. S. Simple Models Show the General Advantages of Dendrites in Coincidence Detection. Journal of Neurophysiology 2007, 97 (5), 3449–3459. https://doi.org/10/dc8tr4.
"""
const MOUSE = Physiology(200Ω * cm, 1700Ω * cm^2, 1μF / cm^2)

## Set synaptic parameters
const human_synapses = TripodSynapses(
    Esyn_soma = Synapse(DuarteGluSoma, MilesGabaSoma),
    Esyn_dend = Synapse(EyalGluDend, MilesGabaDend),
)

const ampa_only = TripodSynapses(
    Esyn_soma = Synapse(DuarteGluSoma, MilesGabaSoma),
    Esyn_dend = Synapse(DuarteGluSoma, MilesGabaDend),
)

const ampa_equivalent = TripodSynapses(
    Esyn_soma = Synapse(DuarteGluSoma, MilesGabaSoma),
    Esyn_dend = Synapse(EyalGluDend_AMPA, MilesGabaDend),
)

const mouse_synapses = TripodSynapses(
    Esyn_soma = Synapse(DuarteGluSoma, MilesGabaSoma),
    Esyn_dend = Synapse(DuarteGluDend, MilesGabaDend),
)

## Set dendritic geometry
const H_ss = (ds = (0, 0), syn_model = human_synapses, species = "H", soma_only = true)
const H_distal_distal = (ds = (400, 400), syn_model = human_synapses, species = "H")
const H_proximal_proximal = (ds = (150, 150), syn_model = human_synapses, species = "H")
const H_distal_proximal = (ds = (400, 150), syn_model = human_synapses, species = "H")
const M_distal_distal = (ds = (400, 400), syn_model = mouse_synapses, species = "M")
const M_proximal_proximal = (ds = (150, 150), syn_model = mouse_synapses, species = "M")
const M_distal_proximal = (ds = (400, 150), syn_model = mouse_synapses, species = "M")

const default_model = H_distal_proximal
const default_synapses = human_synapses
const models = [H_distal_distal, H_proximal_proximal, H_distal_proximal, H_ss]
const labels = ["distal-distal", "proximal-proximal", "distal-proximal", "soma only"]
