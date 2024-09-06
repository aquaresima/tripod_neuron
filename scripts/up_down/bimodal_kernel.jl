using Statistics
using StatsBase
# Find bimodal value
# Here I use a simple algorithm that is described in :
# Journal of the Royal Statistical Society. Series B (Methodological)
# Using Kernel Density Estimates to Investigate Multimodality
# https://www.jstor.org/stable/2985156

# It consists in  using Normal kernels to approximate the data and then leverages a theorem on decreasing monotonicity of the number of maxima as function of the window span.

# Kernel Density Estimation
function KDE(t::Real, h::Int64, data)
    ndf(x, h) = exp(-x^2 / h)
    1 / length(data) * 1 / h * sum(ndf.(data .- t, h))
end

# Distribution
function globalKDE(h::Int64, data; v_range = collect(-90:-35))
    kde = zeros(Float64, length(v_range))
    @fastmath @inbounds for n = 1:length(v_range)
        kde[n] = KDE(v_range[n], h, data)
    end
    return kde
end

#Get its maxima
function get_maxima(data)
    arg_maxima = []
    for x = 2:length(data)-1
        (data[x] > data[x-1]) && (data[x] > data[x+1]) && (push!(arg_maxima, x))
    end
    return arg_maxima
end

#Trash spurious values (below 30% of the true maximum)
function isbimodal(kernel, ratio)
    maxima = get_maxima(kernel)
    z = maximum(kernel[maxima])
    real = []
    for n in maxima
        m = kernel[n]
        if (abs(m / z) > ratio)
            push!(real, m)
        end
    end
    if length(real) > 1
        return true
    else
        return false
    end
end

#Trash spurious values (below 30% of the true maximum)
function count_maxima(kernel, ratio)
    maxima = get_maxima(kernel)
    z = maximum(kernel[maxima])
    real_maxima = []
    for n in maxima
        m = kernel[n]
        if (abs(m / z) > ratio)
            push!(real_maxima, m)
        end
    end
    return length(real_maxima)
end

# Return the critical window (hence the bimodal factor)
function critical_window(data; ratio = 0.3, max_b = 50, v_range = collect(-90:-35))
    for h = 1:max_b
        kernel = globalKDE(h, data, v_range = v_range)
        bimodal = false
        try
            bimodal = isbimodal(kernel, ratio)
        catch
            bimodal = false
            @error "Bimodal failed"
        end
        if !bimodal
            return h
        end
    end
    return max_b
end

# Return the critical window (hence the bimodal factor)
function all_windows(data, ratio = 0.3; max_b = 50)
    counter = zeros(max_b)
    for h = 1:max_b
        kernel = globalKDE(h, data)
        counter[h] = count_maxima(kernel, ratio)
    end
    return counter
end
