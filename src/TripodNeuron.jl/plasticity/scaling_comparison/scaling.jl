##
using Plots
using ProgressBars
using StatsBase
using Statistics
using RecipesBase
using BlockArrays
using Random
using Printf
Nee = 100
jeemin = 1.0
jeemax = 50

default(grid = false, guidefontsize = 15)
##

mutable struct Particle{T<:Real}
    μ::T
    ξ::T
    Particle{T}(x) where {T<:Real} = new{T}(convert(T, x), rand(T) - 0.5)
    Particle{T}() where {T<:Real} = new{T}(convert(T, 1), rand(T) - 0.5)
end


Array{Particle}(T::Type, x::Int, y::Int) =
    hcat([[Particle{T}() for _ = 1:x] for _ = 1:y]...)

Particles(T, x::Int, y::Int) = Array{Particle}(T, x, y)
Particles(x::Int, y::Int) = Array{Particle}(Float32, x, y)
function Particles(m::Matrix)
    p = Particles(typeof(m[1]), size(m)...)
    for x in eachindex(m)
        p[x].μ = m[x]
    end
    return p
end

import Base.sum
sum(p::Array{Particle,2}, kwargs...) = sum(getfield.(p, :μ), kwargs...)

Base.:+(p1::Particle, p2::Particle) = p1.μ + p2.μ
Base.:+(p1::Particle, p2::Real) = p1.μ + p2
Base.:+(p2::Real, p1::Particle) = p1.μ + p2

function CRW!(x::Particle; α::Float32, σ::Float32, θ::Float32, x0::Float32)
    x.ξ = (1 - α) * x.ξ + (α) * (rand() - 0.5)
    x.μ += 0.1 * (θ * (x0 - x.μ) + σ * x.ξ)
end

function CRW!(
    x::Array{Particle,2};
    α::Float32,
    σ::Float32,
    θ::Float32,
    bias::Matrix{Float32},
)
    for i = 1:size(x, 1)
        for j = 1:size(x, 2)
            CRW!(x[i, j], α = α, σ = σ, θ = θ, bias = bias[i, j])
        end
    end
    return x
end

## Scaling and simulation

function run_noscaling(w0, bias, steps)
    w = Particles(w0)
    for _ = 1:steps
        for x in eachindex(w)
            update!(w[x], bias[x])
        end
    end
    return getfield.(w, :μ)
end

function scaling_add!(w, sum_w0)
    Δ = (sum(w, dims = 2) .- sum_w0) / size(w, 2)
    for (neuron_post, w0_row) in zip(eachrow(w), eachrow(Δ))
        for neuron_pre in eachindex(neuron_post)
            neuron_post[neuron_pre].μ -= w0_row[1]
        end
    end
end

function scaling_mult!(w, sum_w0)
    Δ = (sum(w, dims = 2) ./ sum_w0)
    for (neuron_post, w0_row) in zip(eachrow(w), eachrow(Δ))
        for neuron_pre in eachindex(neuron_post)
            neuron_post[neuron_pre].μ /= w0_row[1]
        end
    end
end

function run_scaling(
    w::Matrix{Float32},
    scaling::Function;
    α = 0.5f0,
    θ = 1.0f-2,
    bias = nothing,
    σ = 5.0f0,
    total_steps = 1000,
)

    if isnothing(bias)
        bias = zeros(Float32, size(w)...)
        θ = 0.0f0
    end
    sumw0 = sum(w, dims = 2)[:, 1]
    w = Particles(w)
    update!(w, bias) = CRW!(w, α = α, σ = σ, θ = θ, x0 = bias)
    #
    synapses = zeros(total_steps, size(w)...)
    iterations = ProgressBar(1:total_steps)
    @inbounds @fastmath for tt = 1:total_steps
        @simd for i = 1:size(w, 1)
            for j = 1:size(w, 2)
                synapses[tt, i, j] = w[i, j].μ
                update!(w[i, j], bias[i, j])
                (w[i, j].μ < jeemin) && (w[i, j].μ = jeemin)
                (w[i, j].μ > jeemin) && (w[i, j].μ > jeemax)
            end
        end

        if tt % 50 == 0
            scaling(w, sumw0)
        end
    end
    scaling(w, sumw0)
    return getfield.(w, :μ), synapses
end

function compare_scaling(w0; α = α, bias = nothing, seed = 23, time_step = 1000)
    # Random.seed!(seed)
    w_mul, synapses_mul =
        run_scaling(w0, scaling_mult!, α = α, total_steps = time_step, bias = bias)
    Random.seed!(seed)
    w_add, synapses_add =
        run_scaling(w0, scaling_add!, α = α, total_steps = time_step, bias = bias)
    return w_mul, w_add, synapses_mul, synapses_add
end
