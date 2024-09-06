using JLD2, HDF5, DrWatson

# file = joinpath(@__DIR__,"optimal_kies_rate_50.jld")
# fid = read(h5open(file,"r"))

μmem = -55.0 ## target membrane value
name = "_highres"
name = ""

# file = joinpath(@__DIR__,"optimal_IE_μmem$(μmem)_$name.jld2")
file = joinpath(@__DIR__, "optimal_IE_μmem$(μmem)$name.jld2")
if !isfile(file)
    @error "File $file does not exist, run 3_optimal_kies.jl first"
else
    fid = load(file) |> dict2ntuple
    @info "Loaded balance K_EI from $file"
end
file = joinpath(@__DIR__, "optimal_IE_μmem$(μmem)_soma_only.jld2")
if !isfile(file)
    @error "File $file does not exist, run 10_soma_only.jl for soma first"
else
    fid_soma = load(file) |> dict2ntuple
    @info "Loaded balance K_EI from $file"
end

νs = fid.νs
opt_kies = fid.opt_kies
opt_kies_soma = fid_soma.opt_kies

soma_syn_models =
    (ampa_eq = ampa_equivalent, nmda = nmda_soma, kuhn = ampa_kuhn, ampa = human_synapses)

balance_kie_soma = (
    ampa_eq = opt_kies_soma.AMPA_EQ,
    nmda = opt_kies_soma.NMDA,
    kuhn = opt_kies_soma.KUHN,
    ampa = opt_kies_soma.AMPA,
)

balance_kie_rate = (ampa = opt_kies.dend_AMPA, nmda = opt_kies.dend_NMDA)
balance_kie_gsyn = (ampa = opt_kies.istdp_AMPA, nmda = opt_kies.istdp_NMDA)
min_AMPA = fid.min_AMPA
min_NMDA = fid.min_NMDA
try
    balance_name = fid.name
catch
    balance_name = ""
end

balance_path = @__DIR__

# Balanced condition
# r = 17 = 17
# balance_rates[17] = 4.39523712324829
# L = 31 = 31
# balance_models[L] = 400
# balance_kie_rate.nmda[r, L] = 1.5063886519351648
# L = 6 = 6
# balance_models[L] = 150
# balance_rates[17] = 4.39523712324829
# balance_kie_rate.nmda[r, L] = 1.7871067675001522
