
##

p, q = plot(ylabel = "CRW", xlabel = ""), plot(ylabel = "Correlation", xlabel = "Time ")

for α = 0.1:0.1:0.5
    x = Particle{Float32}(1.0)
    w = zeros(Float32, 10000)
    for n = 2:10000
        CRW!(x; α = 1.0f0, σ = 1.0f0, θ = exp(-1 / Float32(α)), bias = 2.0f0)
        w[n] = x.μ
    end
    plot!(p, w, labels = "")
    plot!(q, autocor(w, 1:1:100), labels = α)
end
pp = plot(p, q, layout = (2, 1))
savefig(pp, joinpath(@__DIR__, "figures", "synaptic_input.pdf"))

## Run example with matrix and no bias

update!(w, bias) = CRW!(w, α = 0.5f0, σ = 1.0f0, θ = 0.0f0, bias = bias)

bias = BlockArray(zeros(Float32, 1, 40), [1], [20, 20])
bias[Block(1, 1)] = fill!(zeros(1, 20), 5)
bias[Block(1, 2)] = fill!(zeros(1, 20), -5)
w0 = rand(1, 40)
w = run_noscaling(w0, bias)
title = round(corspearman(w0[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p1 = scatter(w0', c = :black, markersize = 7, label = "Initial condition")
p1 = scatter!(
    getfield.(w, :μ)',
    label = "Final state",
    legendfontsize = 10,
    markersize = 6,
    title = title,
)

w0 = Matrix(bias)
w = run_noscaling(w0, bias)
title = round(corspearman(w0[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p2 = scatter(w0', c = :black, markersize = 7, label = "")
p2 = scatter!(
    getfield.(w, :μ)',
    label = "",
    legendfontsize = 12,
    markersize = 6,
    xlabel = "Pre-synaptic neurons",
    title = title,
)

bias[1, :] = rand(1, 40)
w0 = zeros(1, 40);
w0[1, :] .= [5sin(x * 2π / 20) for x = 1:40];
w = run_noscaling(w0, bias)
title = round(corspearman(w0[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p3 = scatter(w0', c = :black, markersize = 7, label = "")
p3 = scatter!(
    getfield.(w, :μ)',
    label = "",
    legendfontsize = 12,
    markersize = 6,
    title = title,
)


pp = plot(p1, p2, p3, layout = (1, 3), yaxis = false, yticks = :none, ylims = (-10, 10))

savefig(pp, joinpath(@__DIR__, "figures", "array_no_bias.pdf"))
pp
## Run example with matrix and bias

update!(w, bias) = CRW!(w, α = 0.2f0, σ = 1.0f0, θ = 1.0f0, bias = bias)

bias = BlockArray(zeros(Float32, 1, 40), [1], [20, 20])
bias[Block(1, 1)] = fill!(zeros(1, 20), 5)
bias[Block(1, 2)] = fill!(zeros(1, 20), -5)
w0 = rand(1, 40)
w = run_noscaling(w0, bias)
title = round(corspearman(bias[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p1 = scatter(w0', c = :black, markersize = 7, label = "Initial condition")
p1 = scatter!(getfield.(w, :μ)', label = "Final state", legendfontsize = 10, markersize = 7)
p1 = plot!(bias', c = :blue, markersize = 7, label = "Bias", title = title)

bias[1, :] = [5sin(x * 2π / 20) for x = 1:40]
w = run_noscaling(w0, bias)
title = round(corspearman(bias[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p2 = scatter(w0', c = :black, markersize = 7, label = "")
p2 = scatter!(
    getfield.(w, :μ)',
    label = "",
    legendfontsize = 10,
    markersize = 6,
    xlabel = "Pre-synaptic neurons",
)
p2 = plot!(bias', c = :blue, markersize = 7, label = "", title = title)

bias[1, :] = rand(1, 40)
w0 = zeros(1, 40);
w0[1, :] .= [5sin(x * 2π / 20) for x = 1:40];
w = run_noscaling(w0, bias)
title = round(corspearman(bias[1, :], getfield.(w, :μ)[1, :]), digits = 2)
p3 = scatter(w0', c = :black, markersize = 7, label = "")
p3 = scatter!(getfield.(w, :μ)', label = "", markersize = 7)
p3 = plot!(bias', c = :blue, markersize = 7, label = "", title = title)


pp = plot(p1, p2, p3, layout = (1, 3), yaxis = false, yticks = :none)

savefig(pp, joinpath(@__DIR__, "figures", "array_bias.pdf"))
pp
