Tripod Neuron
-------------

The code presented here simulates a three-compartments neuron model. This package serves for reproducibility and future work based on the journal paper:

__"The Tripod neuron: a minimal structural reduction of the
dendritic tree"__ 
AQ et al.

The TripodNeuron integrates dendritic and soma dynamics with a circuital approximation (passive compartments separated by a resistance). The integration follows the Heun method (second-degree Euler) for each compartment.

## Parameters
The Tripod code collects all the parameters in 4 structs:
  1. < Synapse >
  2. < Receptor >
  2. < AdExParams >  
  2. < PassiveMembraneParameters > 

The parameters are instantiated within the module as constant variable, but can also be defined at runtime.

## Simulation loop
The simulation is run by the function `run_tripod`. 
This function accepts the following parameters:
``` 
# Arguments
  model: dendritic compartment length; ex. (150, 400)
	inputs: Matrix or Function for the inputs, see below
	simtime : duration of the simulation in seconds; ex. 1000

# Keyworded arguments
	species: species for membrane parameters
	ds: dendritic length, overwrites the model argument; ex. (150, 400,)
	synapses= Vector with synaptic efficacies of length six; ex. ones(6) 
	syn_model: Synapse struct with parameters
	AdEx: AdExParams struct with parameters
	postspike: PostSpike struct with parameters
	ext_currents: External currents, length 3. :: zeros(3)
	do_spikes=true: Spike or non-spike somatic model
	soma_only=false: Integrate the model without dendritic compartments

# Recording parameters.
	adapt=false,
	in_curr = false,
	kwargs...
	)
```

When the model is integrated in the `soma_only` configuration, the soma receives synaptic currents from two independent synaptic variables. Which is formally equivalent to dendrites with length zero.

The simulation timestep `dt` is defined as the first variable of the model, several objects depend on it. If you change `dt` consider to restart the Julia REPL. The simulator has not been tested with `dt` different from 0.1 ms

## Inputs:
There are two ways of defining the inputs for the Tripod neuron. The `run_tripod` function will accept this type:
`Union{Matrix{Bool}, Matrix{Int}, Matrix{Float32}, Vector{<:Real}, Function}`

### Matrix inputs
The matrix define all the spikes carried on the Tripod at each moment in the simulation. The dimension of the matrix are fixed. The rows indicate the compartment which is targeted by the spike and the type of neurotransmitters simulated (Glutamatergic/ GABAergic). The column size has to match with the simulated time-steps,(i.e., `simtime/dt`). 

### Function inputs
When the inputs are given in the form of a function, the simulator will call this function for every compartment at each time step.
The function has two mandatory arguments: 
```
function input_function(
      tt::Int64 # The simulation timestep
      compartment::Int64 #The target comapartment.
)

    if tt == 10/dt # at timestep equal to 10 ms
      if compartment == 2 # target the first dendritic compartment.
        return 5. # trigger a spike of synaptic weight equalt to 5 g_syn.
      end
    end

end

```

### Target compartments:
These are the target compartment and the corresponding neurotransmitters.

	# 1-> soma excitatory
	# 2-> d1 excitatory
	# 3-> d2 excitatory
	# 4-> soma inhibitory
	# 5-> d1 inhibitory
	# 6-> d2 inhibitory

## Protocols
The file `simulator/protocols.jl` contains several utilities to generate spike trains and measure the neuron activit. The names of the function are explanatory, they are largely used in the code to reproduce the figures.


## Code structure
The code is structured within these files:

  * _src/equations.jl_ mathematical implementation of model equations
  * _src/simulation.jl_ simulation loop 
  * _src/simulation_sequence.jl_ simulation loop with feedback inhibition

  * _src/simulation_sequence.jl_ simulation loop with feedback inhibition

  * _model/units.jl_: SI units and physiological parameters
  * _model/synapses.jl_ structs and parameters.
  * _model/soma.jl_ AdExp parameters.
  * _model/soma.jl_ AdExp parameters.
  * _model/dendrites.jl_: Struct and constructor for passive dendritic compartments.
