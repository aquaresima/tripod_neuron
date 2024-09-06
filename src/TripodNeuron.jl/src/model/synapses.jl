# Copyright (c) 2022 Alessio Quaresima
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

abstract type AbstractReceptor end

@with_kw struct Receptor <: AbstractReceptor
    E_rev::Float32 = 0.0 # Reversal potential
    τr::Float32 = -1.0f0 # rise timescale
    τd::Float32 = -1.0f0 # decay timescale
    g0::Float32 = 0.0f0 # g0 is peak conductance 
    ## Parameters computed for the double_exp integration
    gsyn::Float32 = g0 > 0 ? g0 * norm_synapse(τr, τd) : 0.0f0
    α::Float32 = α_synapse(τr, τd)
    τr⁻::Float32 = 1 / τr > 0 ? 1 / τr : 0.0f0
    τd⁻::Float32 = 1 / τd > 0 ? 1 / τd : 0.0f0
end

## NMDA Receptor parameters
Mg_mM = 1.0f0
nmda_b = 3.36       #(no unit) parameters for voltage dependence of nmda channels
nmda_k = -0.077     #Eyal 2018

## Voltage dependent receptor (NMDA)
@with_kw struct ReceptorVoltage <: AbstractReceptor
    E_rev::Float32 = 0.0 # Reversal potential
    τr::Float32 = -1.0f0 # rise timescale
    τd::Float32 = -1.0f0 # decay timescale
    g0::Float32 = 0.0f0 # g0 is peak conductance 
    ## Parameters computed for the double_exp integration
    gsyn::Float32 = g0 > 0 ? g0 * norm_synapse(τr, τd) : 0.0f0
    α::Float32 = α_synapse(τr, τd)
    b::Float32 = nmda_b
    k::Float32 = nmda_k
    mg::Float32 = Mg_mM
    τr⁻::Float32 = 1 / τr > 0 ? 1 / τr : 0
    τd⁻::Float32 = 1 / τd > 0 ? 1 / τd : 0
end

struct Synapse
    AMPA::Receptor
    NMDA::ReceptorVoltage
    GABAa::Receptor
    GABAb::Receptor
end

struct Glutamatergic
    AMPA::Receptor
    NMDA::ReceptorVoltage
end

struct GABAergic
    GABAa::Receptor
    GABAb::Receptor
end

@with_kw struct TripodSynapses
    Esyn_dend::Synapse
    Esyn_soma::Synapse
end


function Synapse(glu::Glutamatergic, gaba::GABAergic)
    return Synapse(glu.AMPA, glu.NMDA, gaba.GABAa, gaba.GABAb)
end


#=========================================
			Synaptic fit
=========================================#

function norm_synapse(synapse::Union{Receptor,ReceptorVoltage})
    norm_synapse(synapse.τr, synapse.τd)
end

# Factor used to compute the effective peak conductance. 
# Roth A & van Rossum MCW (2009), Modeling synapses. In Computational Modeling Methods for
# 901 Neuroscientists, ed. De Schutter E, MIT Press, Cambridge, MA, USA, pp. 139–160.
function norm_synapse(τr, τd)
    p = [1, τr, τd]
    t_p = p[2] * p[3] / (p[3] - p[2]) * log(p[3] / p[2])
    return 1 / (-exp(-t_p / p[2]) + exp(-t_p / p[3]))
end

# α is the factor that has to be placed in-front of the differential equation as such the analytical integration corresponds to the double exponential function. Further details are discussed in the Julia notebook about synapses
function α_synapse(τr, τd)
    return (τd - τr) / (τd * τr)
end


#=====================================
#  			Synapses params
=====================================#

"""
Christof Koch. Biophysics of Computation: Information Processing in Single Neurons, by.Trends in Neurosciences, 22842(7):328–329, July 1999. ISSN 0166-2236, 1878-108X. doi: 10.1016/S0166-2236(99)01403-4.
"""
const KochGlu =
    Glutamatergic(Receptor(E_rev = 0.00, τr = 0.2, τd = 25.0, g0 = 0.73), ReceptorVoltage())


"""
Guy Eyal, Matthijs B. Verhoog, Guilherme Testa-Silva, Yair Deitcher, Ruth Benavides-Piccione, Javier DeFelipe, Chris-832tiaan P. J. de Kock, Huibert D. Mansvelder, and Idan Segev. Human Cortical Pyramidal Neurons: From Spines to833Spikes via Models.Frontiers in Cellular Neuroscience, 12, 2018. ISSN 1662-5102. doi: 10.3389/fncel.2018.00181.
"""

const EyalGluDend = Glutamatergic(
    Receptor(E_rev = 0.0, τr = 0.26, τd = 2.0, g0 = 0.73),
    ReceptorVoltage(E_rev = 0.0, τr = 8, τd = 35.0, g0 = 1.31),
)


"""
Richard Miles, Katalin Tóth, Attila I Gulyás, Norbert Hájos, and Tamas F Freund.  Differences between Somatic923and Dendritic Inhibition in the Hippocampus.Neuron, 16(4):815–823, April 1996. ISSN 0896-6273. doi: 10.1016/924S0896-6273(00)80101-4.
"""
const MilesGabaDend = GABAergic(
    Receptor(E_rev = -75.0, τr = 4.8, τd = 29.0, g0 = 0.27),
    Receptor(E_rev = -90.0, τr = 30, τd = 400.0, g0 = 0.006),
)

const MilesGabaSoma =
    GABAergic(Receptor(E_rev = -75.0, τr = 0.3, τd = 15.0, g0 = 0.38), Receptor())


"""
Renato Duarte and Abigail Morrison. Leveraging heterogeneity for neural computation with fading memory in layer 2/3808cortical microcircuits.bioRxiv, December 2017. doi: 10.1101/230821.
"""
const DuarteGluSoma = Glutamatergic(
    Receptor(E_rev = 0.0, τr = 0.25, τd = 2.0, g0 = 0.73),
    ReceptorVoltage(E_rev = 0.0),
)

const DuarteGluDend = Glutamatergic(
    Receptor(E_rev = 0.0, τr = 0.25, τd = 2.0, g0 = 0.73),
    ReceptorVoltage(E_rev = 0.0, τr = 0.99, τd = 100.0, g0 = 0.159),
)

const DuarteGabaSoma = GABAergic(
    Receptor(E_rev = -75.0, τr = 0.5, τd = 6.0, g0 = 0.265),
    Receptor(E_rev = -90.0, τr = 30, τd = 100.0, g0 = 0.006),
)

## AMPA equivalent
## The AMPA g0 is set to 2.0 such that the EPSP is comparable
const EyalGluDend_AMPA =
    Glutamatergic(Receptor(E_rev = 0.0, τr = 0.26, τd = 2.0, g0 = 2.0), ReceptorVoltage())
