module TripodNeuron
using Random
using Distributions
using Logging
using Printf
using Unitful
using Parameters
using UnPack

include("model/parameters.jl")
include("simulator/equations.jl")
include("simulator/simulation.jl")
include("simulator/simulation_feedback.jl")
include("simulator/protocols.jl")

## updown:
include("simulator/balance.jl")
include("simulator/vogels_istdp.jl")
include("../balance_EI/load.jl")

TN = TripodNeuron
export TN

end
