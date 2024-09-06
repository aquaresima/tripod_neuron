using RollingFunctions
using StatsBase

##6. Stern, E. A., Kincaid, A. E. & Wilson, C. J. Spontaneous subthreshold membrane
#To measure the duration of each state transition, two thresholds were set at one fourth and three fourths of the distance between the peaks of the membrane potential distribution.

"""
function voltage_cross_counter(v, threshold_up, threshold_down)
	For each value of v, verifies if one of the threshold has been crossed. Create an array with the length of the interval between each crossing. The array is assigned to the up or down state depending on the threshold crossed.
"""
function voltage_cross_counter(v, threshold_up = 3 / 4, threshold_down = 1 / 4)
    v = cap_voltage(v)
    _v = v[10_000:end-1]
    threshold_up = threshold_up * (maximum(_v) - minimum(_v)) + minimum(_v)
    threshold_down = threshold_down * (maximum(_v) - minimum(_v)) + minimum(_v)
    v = cap_voltage(v)
    @info "Threshold for up-down states are" threshold_down, threshold_up
    up_states = Vector{Vector{Int64}}(undef, 0)
    down_states = Vector{Vector{Int64}}(undef, 0)
    time_points = Vector{Int64}(undef, 0)
    state = :none
    for i in eachindex(v)
        if v[i] >= threshold_up
            if state == :down
                # @warn v[i]
                # @warn "Down state is too short: $(length(time_points))"
                push!(down_states, time_points)
                @assert all(v[time_points] .< threshold_up)
                time_points = Vector{Int64}(undef, 0)
            end
            state = :up
        elseif v[i] <= threshold_down
            # @show time_points
            if state == :up
                # @warn "Up state is too short: $(length(time_points))"
                push!(up_states, time_points)
                # @assert all(v[time_points] .> threshold_down)
                time_points = Vector{Int64}(undef, 0)
            end
            state = :down
        end
        @assert (state == :up || state == :down || state == :none)
        if state == :up || state == :down
            push!(time_points, i)
        end
        # if state == :up
        # 	@assert v[i] > threshold_down
        # elseif state == :down
        # 	@assert v[i] < threshold_up
        # end
    end
    return (up_states, down_states, threshold_down, threshold_up)
    # return v, threshold_down, threshold_up
end


"""
Measure the membrane potential average, the duration in milli seconds, the number of spikes, and standard deviation for each up and down state.
Return:
	- up: a matrix with the average, duration, standard deviation, and number of spikes for each up state
	- down: a matrix with the average, duration, standard deviation, and number of spikes for each down state

Order of the columns:
	- duration
	- average
	- standard deviation
	- number of spikes

"""
function updown_analysis(voltage)
    up, down, _, _ = voltage_cross_counter(voltage)
    av_up = Vector{}()
    av_down = Vector{}()
    voltage = Float64.(voltage)
    for i in eachindex(up)
        n_spikes = round(Int, TN.get_spike_rate(voltage[up[i]]) * 10000 / length(up[i]))
        # @warn "Up spikes" n_spikes
        xx = up[i]
        vv = cap_voltage(voltage[xx])
        push!(av_up, [length(xx) * 0.1, mean(vv), std(vv), n_spikes])
    end
    for i in eachindex(down)
        n_spikes = round(Int, TN.get_spike_rate(voltage[down[i]]) * 10000 / length(down[i]))
        xx = down[i]
        vv = cap_voltage(voltage[xx])
        push!(av_down, [length(xx) * 0.1, mean(vv), std(vv), n_spikes])
    end
    # return up,down
    return (up = hcat(av_up...), down = hcat(av_down...))

end

function cum_prob(x)
    """
    Compute the cumulative probability of x
    """
    xx = sort(x)
    xx = xx ./ sum(xx)
    return cumsum(xx)
end

function cum_prob(x)
    return cumsum(x ./ sum(x))
end
