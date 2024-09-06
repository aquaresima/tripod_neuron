using DrWatson, Logging
@quickactivate "Tripod"
using TripodNeuron
using Statistics

function get_input_cv(β, rate, τ)
    cond = TN.get_balance_conditions(1, 1, rate)
    spikes = TN.make_spikes(10_000, β; τ = τ, rate = cond.rate)[3, 20_000:end]
    isi = Vector{Float64}()
    last = 0
    for x = 1:length(spikes)
        ns = spikes[x]
        z = x
        while ns > 0
            ns -= 1
            push!(isi, z * TN.dt - last)
            last = z * TN.dt
            z += 1
        end
    end
    return TN.CV_isi(Float32.(isi))
end



rates = 1:5:length(TN.νs)-5
rates = 1:5:30
βs = 0:20:1000
τs = collect(exp.(range(log(0.1), log(1000), length = 50)))
cvs = zeros(length(τs), length(βs), length(rates))
for r in eachindex(rates[1:1])
    for n in eachindex(τs)
        @info "τ: $(τs[n]), rate: $(rates[r])"
        Threads.@threads for m in eachindex(βs)
            cvs[n, m, r] = mean(get_input_cv(βs[m], rates[r], τs[n]) for x = 1:50)
        end
    end
end

data = (cvs = cvs, rates = rates, βs = βs, τs = τs)
safesave(datadir("up_down", "inputs", "cv_inputs_norate.bson"), @dict data)
