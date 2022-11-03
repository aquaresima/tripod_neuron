using RollingFunctions
using StatsBase

@unpack AP_membrane = PostSpike()
@show AP_membrane
function get_spikes(voltage::Array{Float32,3}; column=1)
    """
    Voltage: (cells, compartments, time)
    column : 1 or 3
    """
    if !(column ∈ [1,3])
        throw(DomainError("column must be 1 (cells to columns) or 3 (time to columns)"))
    end
    if column == 1
        spikes = falses(size(voltage,1),size(voltage,3))
        for n in 1:size(voltage,1)
            spikes[n,:] .= [ v ==AP_membrane for v in voltage[n,1,:]]
        end
    end
    if column == 3
        spikes = falses(size(voltage,3),size(voltage,1))
        for n in 1:size(voltage,1)
            spikes[:, n] .= [ v == AP_membrane for v in voltage[n,1,:]]
        end
    end
    return spikes
end

get_spikes(voltage::Array{Float32,1})::Array{Bool,1} = [v == AP_membrane for v in voltage[:]]

function get_spikes(voltage::Array{Float32,2})
    spikes = falses(size(voltage,2),size(voltage,1))
    for n in 1:size(voltage,1)
        spikes[:,n] .= [ v == AP_membrane for v in voltage[n,:]]
    end
    return spikes
end

##
get_spike_times(voltage::Array{Float32,2}) = findall(x->x==AP_membrane, get_spikes(voltage))

function get_spike_times(spikes::BitArray{2})
    times= Vector()
    for neuron in eachrow(spikes)
        push!(times,findall(neuron))
    end
    return times
end

function get_spike_times(spikes::Matrix{Int})
    times= Vector()
    for neuron in eachrow(spikes)
        push!(times,findall(x->x>0,neuron))
    end
    return times
end
##
get_spike_times(spikes::BitArray{1}) = findall(spikes)
get_spike_times(voltage::Array{Float32,1}) = findall(get_spikes(voltage))

##
get_spike_rate(v::Array{Float32,2})::Float32 = float(sum(v[1,:].== AP_membrane)/length(v[1,:])*10000)
get_spike_rate(v::Array{Float32,1}) = sum(v.== AP_membrane)/length(v)*10000
get_spike_rate(v::Array{Float32,3}) = [sum(v[i,1,:].== AP_membrane)/length(v[i,1,:])*10000 for i in 1:size(v,1)]
get_spike_rate(spikes::BitArray{1}) = sum(spikes)/length(spikes)*10000
get_spike_rate(spikes::BitArray{2}) = map(x->get_spike_times(x), eachrow(spikes))

##
function get_isi(input::Union{Matrix{Float32}, Vector{Float32}, BitArray{1}})
    if isa(input,Array{Float32,2})
        return diff(get_spike_times(input[1,:]))
    elseif isa(input,BitArray{1})
        return diff(get_spike_times(input))
    elseif isa(input,Array{Float32,1})
        return diff(get_spike_times(input))
    end
end

function CV_isi(intervals::Vector{Float32})
    # input::Union{Matrix{Float32}, Vector{Float32}, BitArray{1}})
    _cv = sqrt(var(intervals)/mean(intervals)^2)
    return isnan(_cv) ? 0. : _cv
end

function CV_voltage(voltage::Vector{Float32})
    intervals = get_isi(voltage)
    _cv = sqrt(var(intervals)/mean(intervals)^2)
    return isnan(_cv) ? 0. : _cv
end

# get_cv(v) =
function get_cv(;v=nothing, isi=nothing)
    if !isnothing(v)
        return CV_voltage(v)
    elseif !isnothing(isi)
        return CV_isi(isi)
    end
end


"""
Function to estimate the istantaneuous firing rate of a neuron.
Calculate the number of spikes every millisecond and make bins, then scale up to the firing rate in second

Parameters
==========

"""
function firing_rate(spikes::BitArray{2}; time_step=dt, window=10 )
    cells, duration = size(spikes)
    scale = 1000 / window / time_step
    window = round(Int,window /time_step)
    Δt = window
    rate = zeros(cells, round(Int,duration/Δt))
    gaussian(x::Real) = exp(-(x-window/2)^2/2window)
    weights=gaussian.(0:9)
    # weights[1:round(Int, Δt/2)].=0
    weights /= sum(weights)
    for (n,cell) in enumerate(eachrow(spikes))
        for x in 1:size(rate)[2]-1
            # rate[n, x]= mean(running(sum, spikes[n,:],window, weights)[1+Δt*(x-1):Δt*x])
            rate[n, x]= sum(spikes[n,1+Δt*(x-1):Δt*x])
        end
        rate[n, :] .= runmean(rate[n,:],round(Int,1/time_step), weights)
    end
    return rate * scale
end


function firing_rate(spikes::BitArray{1}; time_step=dt, window=10 )
    my = falses(1, size(spikes)[1])
    my[1,:] .= spikes
    return firing_rate(spikes; time_step=time_step, window=window )[1,:]
end

# function firing_rate(spikes::BitArray{2}; time_step=dt, slide_offset::Int64 = 20, window::Int64 = 50)
#     ### the time of each bin is 0.1 ms.
#     cells, duration = size(spikes)
#     max_duration(x::Int)= x > duration ? duration : x
#     ## make integer values
#     slide_offset = ceil(Int,slide_offset/time_step)
#     window  = round(Int,window/time_step)
#     time_points = floor(Int,duration/slide_offset)
#     rate = zeros(cells, time_points)
#     for n in 1:time_points
#         x = 1+(n-1)*slide_offset
#         # slide*(slide_offset+1):slide:slide*time_points)
#         rate[:, n]=sum(spikes[:,x:max_duration(x+window)],dims=2)*(1000/window)
#         # *(1000/coarse)
#     end
#     return rate
# end



"""
Average membrane potential

Parameters
==========

`voltage`: [neurons, time] matrix with membrane voltage
`window`: length of the averaging window, in ms
`slide_offset`: shift the window of ms for each new point
"""

function average_potential(membrane::Array{Float32,2}, time_step::Float32 ; slide_offset::Float32 = 20, window::Float32 = 50)
    ### the time of each bin is 0.1 ms.
    cells, duration = size(membrane)
    max_duration(x::Int)= x > duration ? duration : x
    ## make integer values
    slide_offset = ceil(Int,slide_offset/time_step)
    window  = round(Int,window/time_step)
    time_points = floor(Int,duration/slide_offset)
    average = zeros(time_points, cells)
    for n in 1:time_points
        x = 1+(n-1)*slide_offset
        # slide*(slide_offset+1):slide:slide*time_points)
        average[n, :]=mean(membrane[:,max_duration.(x:round(Int,window/10):x+window)],dims=2)
    end
    return average
end
