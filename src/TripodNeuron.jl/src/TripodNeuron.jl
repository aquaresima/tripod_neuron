module TripodNeuron
using Random
using Distributions
using Logging
using Printf
using Unitful
using Parameters
using UnPack

include("model/parameters.jl")
include("equations.jl")
include("simulator/simulation.jl")
include("simulator/simulation_feedback.jl")
include("simulator/protocols.jl")


TN = TripodNeuron
export TN

end
