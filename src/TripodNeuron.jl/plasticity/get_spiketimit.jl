
## Load SpikeTimit

path_spiketimit = "/home/cocconat/Documents/Research/phd_project/speech/Spike TIMIT/src/SpikeTimit.jl"
include(path_spiketimit)
path = "/home/cocconat/Documents/Research/phd_project/speech/litwin-kumar_model_thesis/Spike TIMIT"
path_data = "/home/cocconat/Documents/Research/phd_project/speech/datasets/Spike TIMIT"
#Create the path strings leading to folders in the data set
train_path = joinpath(path_data, "train");
dict_path = joinpath(path_data, "DOC/TIMITDIC.TXT");
train = SpikeTimit.create_dataset(; dir = train_path)
dict = SpikeTimit.create_dictionary(file = dict_path)
