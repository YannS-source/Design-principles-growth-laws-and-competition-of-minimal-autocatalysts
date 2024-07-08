import numpy as np
import pickle

path_where_data_is = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/Data/Figure2"

####

with open(path_where_data_is + "/all_epsilons", "rb") as fp:   # Unpickling
    all_epsilon = pickle.load(fp)
with open(path_where_data_is + "/mean_efficiency_5", "rb") as fp:   # Unpickling
    mean_efficiency_for_a_we_barrier_5 = pickle.load(fp)
with open(path_where_data_is + "/error_efficiency_5", "rb") as fp:   # Unpickling
    error_efficiency_for_a_we_barrier_5 = pickle.load(fp)
with open(path_where_data_is + "/mean_efficiency_10", "rb") as fp:   # Unpickling
    mean_efficiency_for_a_we = pickle.load(fp)
with open(path_where_data_is + "/error_efficiency_10", "rb") as fp:   # Unpickling
    error_efficiency_for_a_we = pickle.load(fp)

with open(path_where_data_is + '/saved_dictionary.pkl', 'rb') as f:
    dic_are_in = pickle.load(f)

###

all_barriers = [5, 10, 15, 20]
theory_epsilons = np.arange(0, 22, 0.1)


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)
plt.xticks(theory_epsilons , theory_epsilons)

# plot the efficiencies obtained from the MD simulations
plt.errorbar(all_epsilon, mean_efficiency_for_a_we_barrier_5, yerr=error_efficiency_for_a_we_barrier_5, capsize=10,capthick=3,fmt='o', color='#2187bb')
plt.errorbar(all_epsilon, mean_efficiency_for_a_we, yerr=error_efficiency_for_a_we, fmt='o', capsize=10,capthick=3, color = "#ff7f0e")

# plot the value according to Table S1, Comprehensive Markov model
plt.plot(dic_are_in[0].keys(), dic_are_in[0].values(),  linewidth=7, color="#2187bb", label=r'$ \epsilon_{AB}=$' + r"${}$".format(all_barriers[0]))
plt.plot(dic_are_in[0].keys(), dic_are_in[1].values(),  linewidth=7, color="#ff7f0e",label=r'$ \epsilon_{AB}=$' + r"${}$".format(all_barriers[1]))
plt.plot(dic_are_in[1].keys(), dic_are_in[2].values(),  linewidth=7, color="#2ca02c", label=r'$ \epsilon_{AB}=$' + r"${}$".format(all_barriers[2]))
plt.plot(dic_are_in[2].keys(), dic_are_in[3].values(),  linewidth=7, color="#9467bd", label=r'$ \epsilon_{AB}=$' + r"${}$".format(all_barriers[3]))
plt.xticks(size = 35)
plt.yticks(size = 35)
labels = [item.get_text() for item in ax.get_xticklabels()]
ax.set_xticklabels([str(int(round(float(label), 0))) for label in labels])
plt.locator_params(axis='x', nbins=6)
plt.locator_params(axis='y', nbins=5)
# Get locations and labels
plt.ylim(ymin=None, ymax=None)
locs, labels = plt.yticks()
plt.legend(prop={'size':35})
plt.show()


