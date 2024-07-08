import numpy as np
import ast

################

"""
Create the df from the txt files of the run, and plot
"""
path_serv = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS1/Increase_Barrier"
all_run = np.arange(0, 201, 1)
all_energetic_barriers = np.arange(0, 12, 1)
all_epsilon =  [0]
all_beta = [0.15]
a_beta = all_beta[0]
all_rcut = [1.1]
diameter_particle = 1


list_all_run = []
list_all_epsilon = []
list_all_tilts = []
list_all_release_time = []
list_all_rcut = []
list_all_type = []
list_all_energetic_barriers = []

type_of_reaction = ['alpha_L4']
for a_run in all_run:
    for a_type in type_of_reaction:
        for an_energetic in all_energetic_barriers:
            if a_type == 'alpha_L4':
                list_all_type.append("alpha")
                try:
                    file = path_serv + "/alpha_L4" + "/" + 'length=4_diameterMain=1_r_cut=1.1_Ebarrrier={}_run={}.txt'.format(an_energetic, a_run)
                    with open(file,
                              'r') as f:
                        list_release = ast.literal_eval(f.read())
                    list_all_energetic_barriers.append(an_energetic)
                    list_all_run.append(a_run)
                    list_all_release_time.append(list_release[2])
                    list_all_type.append("alpha")
                except FileNotFoundError:
                    print("Wrong file or file path", a_type, a_run)
                    print(file)
                    print("/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Energetic_old/Model1/Compute_alpha_beta_EnB_v2/alpha/length=4_diameterMain=1_r_cut=1.1_Ebarrrier=0_run=0.txt")




import pandas as pd
df1 = pd.DataFrame(list(zip(list_all_type, list_all_energetic_barriers, list_all_run, list_all_release_time)),
               columns =["Type", "Energy Barrier", "Run", "Reaction Time"])



df1_avr = df1.groupby(["Energy Barrier"])["Reaction Time"].mean()
print("The spontaneous reaction, in L = 4:", df1_avr)


v6 = [64.092596*(np.exp((i-6)/(1.1))) for i in all_energetic_barriers]



import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)
import matplotlib.pyplot as plt
import statistics
from math import sqrt
def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.25):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / sqrt(len(values))

    left = x - horizontal_line_width / 2
    top = mean - confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean + confidence_interval
    plt.plot([x, x], [top, bottom], color=color)
    plt.plot([left, right], [top, top], color=color)
    plt.plot([left, right], [bottom, bottom], color=color)
    plt.plot(x, mean, 'o', color='#1f77b4')

    return mean, confidence_interval


plt.xticks(all_energetic_barriers, all_energetic_barriers)
plt.plot(all_energetic_barriers, v6, linewidth = 5, color='#FFA500', label=r'$T_{sp}(h_{AB}^- = 6)e^{h_{AB}^- - 6}$')
for a_barrier in all_energetic_barriers:
    dd = df1.loc[df1["Energy Barrier"] == a_barrier]
    plot_confidence_interval( a_barrier, dd["Reaction Time"], horizontal_line_width=0.5)

plt.xlabel(r"$\varepsilon_{AB}^+$", size=37)
plt.ylabel(r"$T_{A+B \to AB}(L^2=16)$", size=37)
plt.yscale("log")
plt.xticks(size = 35)
plt.yticks(size = 35)
plt.locator_params(axis='x', nbins=5)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='#1f77b4', label='Label1')
black_patch = mpatches.Patch(color='#FFA500', label='Label2')
plt.legend(handles=[red_patch, black_patch], labels=['Simulations', r'$T_{A+B \to AB}(\varepsilon_{AB}^+ = 6)e^{\varepsilon_{AB}^+ - 6}$'], prop={'size':30})
plt.show()









