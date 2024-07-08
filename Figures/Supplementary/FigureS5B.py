import numpy as np
import ast


######################################
######## Spontaneous Reaction ########
######################################

path_serv = '/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS5/Spontaneous'
all_energetic_barriers = np.arange(0, 9, 1)
all_run = np.arange(0, 201, 1)
all_we = [40]
diameter_particle = 1
all_tilt = [0.0]


list_all_run = []
list_all_epsilon = []
list_all_tilts = []
list_all_release_time = []
list_all_rcut = []
list_all_type = []
list_all_energetic_barriers = []
list_all_run_beta = []
list_all_release_time_beta = []

type_of_reaction = ['alpha', 'beta']
for a_run in all_run:
    for a_type in type_of_reaction:
        for an_energetic in all_energetic_barriers:
            if a_type == 'alpha':
                try: # length=4_diameterMain=1_r_cut=1.1_Ebarrrier=7_run=187
                    file = path_serv + "/alpha" + "/" + 'length=4_diameterMain=1_r_cut=1.1_Ebarrrier={}_run={}.txt'.format(an_energetic, a_run)
                    # print(file)
                    # result_run = [a_run, time]
                    with open(file,
                              'r') as f:
                        list_release = ast.literal_eval(f.read())
                    list_all_energetic_barriers.append(an_energetic)
                    list_all_run.append(a_run)
                    list_all_release_time.append(list_release[2])
                    list_all_type.append("alpha")
                except FileNotFoundError:
                    print("Wrong file or file path", an_energetic, a_type, a_run)
                    # print(file)
                    # print("/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Model2_v3/Spontaneous_L4_Model2_v3/alpha/length=4_diameterMain=1_r_cut=1.1_Ebarrrier=0_run=0.txt")
            elif a_type == 'beta':
                list_all_type.append("beta")
                try:
                    file = path_serv + "/beta" + "/" + 'length=4_diameterMain=1_r_cut=1.1_run={}.txt'.format(a_run)
                    with open(file,
                              'r') as f:
                        list_release = ast.literal_eval(f.read())
                    list_all_run_beta.append(a_run)
                    list_all_release_time_beta.append(list_release[1])
                except FileNotFoundError:
                    print("Wrong file or file path", an_energetic, a_type, a_run)



print("MEAN TIME BETA: ", np.mean(list_all_release_time_beta))

import pandas as pd
dfspontaneous = pd.DataFrame(list(zip(list_all_type, list_all_energetic_barriers, list_all_run, list_all_release_time)),
               columns =["Type", "Energy Barrier", "Run", "Reaction Time"])



dfspontaneous_avr = dfspontaneous.groupby(["Energy Barrier"])["Reaction Time"].mean()
print("The spontaneous reaction, in L = 4:", dfspontaneous_avr)



#################
### On Catalyst ##
#################
path_serv = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS5/On_Catalyst"
all_run = np.arange(0, 221, 1)
all_we = [40]
all_energetic_barriers = np.arange(0, 13, 1)
diameter_particle = 1
all_tilt = [0.0]

list_all_run = []
list_all_epsilon = []
list_all_tilts = []
list_all_rcut = []
list_all_energetic_barriers = []
list_spontaneous = []
list_moment_in_3 = []
list_mean_time_in_2 = []
list_time_2_to_3 = []
pb_from_1_to_2 = []
list_first_time_in_2 = []

for a_we in all_we:
    for a_run in all_run:
        for an_energetic in all_energetic_barriers:
            for a_tilt in all_tilt: # the tilt represent the angle for the patches.
                try:
                    file = path_serv + "/" + 'length=4_tilt={}_wepsilon={}_repbarrier={}_run{}_diameterMain=1_2to3.txt'.format(a_tilt, a_we, an_energetic, a_run)
                    # result_run = [a_length, a_tilt, a_weak_epsilon, a_repulsive_barrier, a_run, reaction_time]
                    with open(file,
                              'r') as f:
                        list_release = ast.literal_eval(f.read())

                        list_moment_in_3.append(list_release[5])

                        list_all_epsilon.append(a_we)
                        list_all_energetic_barriers.append(an_energetic)
                        list_all_run.append(a_run)
                        list_all_tilts.append(a_tilt)
                except FileNotFoundError:
                    print("Wrong file or file path", an_energetic, a_run)
                    print(file)



import pandas as pd
df = pd.DataFrame(list(zip(list_all_epsilon, list_all_energetic_barriers, list_all_tilts, list_all_run, list_moment_in_3)),
               columns =["Weak Epsilon", "Energy Barrier", "tilt", "Run", "Moment in 3"])



df1 = df.groupby(["Energy Barrier"])["Moment in 3"].mean()




######################################
######## Plots ######################
######################################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)
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
v10 = [102.957362 * np.exp((i-8)/1) for i in all_energetic_barriers]
exp_curve = [2322.852114 * np.exp(i-6/1) for i in all_energetic_barriers]
plt.plot(all_energetic_barriers, exp_curve, linewidth = 5, color="darkblue", label=r'$T_{A+B \to AB}$($L^2 = 16$)')
plt.plot(all_energetic_barriers, v10, linewidth = 5, color="red", label=r"$T_{AB(A+B) \to AB(AB)}$")
for a_energy in all_energetic_barriers:
    dd = df.loc[df["Energy Barrier"] == a_energy]
    plot_confidence_interval(a_energy, dd["Moment in 3"], horizontal_line_width=0.25)
    dd = dfspontaneous.loc[dfspontaneous["Energy Barrier"] == a_energy]
    if dd.empty == False:
        plot_confidence_interval(a_energy, dd["Reaction Time"], horizontal_line_width=0.25)

print(exp_curve)
plt.xlabel(r'$\epsilon_{AB}^+$', size=37)
plt.ylabel(r'Time', size=37)
plt.locator_params(axis='x', nbins=7)
plt.yscale("log")
plt.xticks(size = 35)
plt.yticks(size = 35)
plt.legend(prop={'size':30})
plt.show()
