import numpy as np
import pickle


all_barriers = np.arange(0, 36, 1)

############
path_where_data_is = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/Data/Figure2"

# plot the value according to Table S1, Comprehensive Markov model
with open(path_where_data_is + "/max_efficiency_first", "rb") as fp:   # Unpickling
    efficiency_max_first = pickle.load(fp)
with open(path_where_data_is + "/max_efficiency_second", "rb") as fp:   # Unpickling
    efficiency_max_second = pickle.load(fp)
with open(path_where_data_is + "/max_efficiency_third", "rb") as fp:   # Unpickling
    efficiency_max_third = pickle.load(fp)

############
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)

plt.plot(all_barriers, efficiency_max_first,  linewidth=7, color = '#8c564b', label=r'$ L^2= 144$')
plt.plot(all_barriers, efficiency_max_second, linewidth=7, color = '#e377c2', label=r'$ L^2= 196$')
plt.plot(all_barriers, efficiency_max_third,  linewidth=7, color = '#7f7f7f', label=r'$ L^2= 256$')

plt.locator_params(axis='x', nbins=6)
plt.locator_params(axis='y', nbins=5)
plt.xticks(size = 35)
plt.yticks(size = 35)
plt.legend(prop={'size':35})
plt.show()