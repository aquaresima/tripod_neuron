"""
Compute inhibitory balance for the tripod model.
"""

include("inhibitory_balance.jl")

nmda = false
μmem = -50.0
@info "Compute balance, NMDA: $nmda, \t μmem = $μmem"
data, νs, ds, kies = compute_balance_condition(μmem, nmda)
data = @strdict data νs ds kies
syn = nmda ? "NMDA" : "AMPA"
file = datadir("up_down", "balance", "inhibitory_balance_high_res_$(syn)_μmem$(μmem).jld2")
mkpath(dirname(file))
safesave(file, data)

μmem = -55.0
@info "Compute balance, NMDA: $nmda, \t μmem = $μmem"
data, νs, ds, kies = compute_balance_condition(μmem, nmda)
data = @strdict data νs ds kies
syn = nmda ? "NMDA" : "AMPA"
file = datadir("up_down", "balance", "inhibitory_balance_high_res_$(syn)_μmem$(μmem).jld2")
mkpath(dirname(file))
safesave(file, data)
