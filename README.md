# Tripod Neuron
-------------

The code presented here simulates a three-compartments neuron model. This package serves for reproducibility and future work based on the journal article (under review):

__"The Tripod neuron: a minimal structural reduction of the
dendritic tree"__


The package is composed of five folders:

- `src`: contains the source code of the Tripod model, with detailed Readme of the package.
- `papers`: contains the files for reproducing all the figures of the paper
- `scripts`: contains utilities
- `data`
- `plots`

In the `papers/dendritic_memory` there are eight folders, one for each figure of the paper. Within it, there are the scripts that are used to reproduce the results and the panels.

Fig6 to Fig9 are obtained with data from relatively long simulations (between 5 mins to 1 hour on server for calculus with 24-threads). In order to facilitate testing of the package, we added pre-computed data in the folder `data/dendritic_memory`, the plot file will load the data from there. 
If you intend to re-compute it, the `DrWatson.safesave` function will automatically swap the files, the newer one will become accessible and the old data file will be renamed with `filename_#1`.

In `scripts` there is one file that is used to generate the style for the plots and one file that automatically run all the scripts in the paper folder.
The `run_all.jl` file include all the files that have plot functions and generates the plots in `plots/dendritic_memory`; it does not run the files marked with `<filename>_run.jl`. These files may require long computation and should be run only if intended.
If the `run_all.jl` script won't work for some unpredicted reason - we suggest to run it from within a REPL-, it is possible to reproduce the figures running each script separatedly. The scripts with the `_long.jl` will take time to run, but can be skipped if the data are in the folder.

The figures generated may present small variations to the ones presented in the paper becasue during the refactoring of the code small variations were made; for example, for consistency across the several stimuli patterns used in the experiments. The results are in agreement with the conclusion drawn in the manuscript.

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
> Tripod

The Julia version used: `julia = "1.8.2"`

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

