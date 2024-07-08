import numpy as np
import ast

################

"""
Create the df from the txt files of the run, and plot
"""

# Data for spontaneous reaction as a function of the area:
path_serv_move = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/Data/Figure2/Spontaneous_Reaction_Move"

all_run = np.arange(0, 221, 1)
energetic_barrier = 0
all_length = np.arange(4, 12, 1)

##################
##################
list_all_run = []
list_all_length = []
list_all_times = []

for a_run in all_run:
    for a_length in all_length:
        try:
            file = path_serv_move + '/length={}_diameterMain=1_r_cut=1.1_run={}.txt'.format(a_length, a_run)
            with open(file,
                      'r') as f:
                list_release = ast.literal_eval(f.read())
            list_all_run.append(a_run)
            list_all_times.append(list_release[2])
            list_all_length.append(a_length)
        except FileNotFoundError:
            print("Wrong file or file path", a_length, a_run)
            print(file)
            print("/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Model1/Tcol_with_inert_catalyst/length=4_wepsilon=0_repbarrier=0_run0_diameterMain=1.txt")

import pandas as pd
df1 = pd.DataFrame(list(zip(list_all_length, list_all_run, list_all_times)),
               columns =["Lengths", "Run", "Reaction Time"])
df1.to_csv("/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Model1/Tcol_increasing_length_inert.csv", index=False)

df1_avr = df1.groupby(["Lengths"])["Reaction Time"].mean()
print("The spontaneous reaction with different L, move:", df1_avr)


################## Theoretical expectation
def T_col_avr(radius_disk=10, rho1=2.35, rho2=2.35, diffusion2d =0.037):
    # the max rho1 + rho2 should be such that rho1 + 2 * rho2 = radius_disk
    Rbar = radius_disk
    p = rho1 + rho2
    if p/3 >= radius_disk:
        print("Particle bigger than volume")
    D = diffusion2d
    pre_factor = 1 / (Rbar -p)
    first_int = (Rbar**2)/(2*D) * ((Rbar)*np.log(Rbar/p) - Rbar + p)
    second_int = (Rbar**3)/(12*D) - (Rbar*p**2)/(4*D) + (2*p**3)/(12*D)
    return pre_factor * (first_int - second_int)
list_reaction_time_theory = []
list_reaction_time_theory_nmove = []
for a_length in all_length:
    radius_from_length = np.sqrt((a_length ** 2) / np.pi)
    a_Tcol = T_col_avr(radius_disk=radius_from_length, rho1=0.6, rho2=0.6,
                       diffusion2d=0.2)
    list_reaction_time_theory.append(a_Tcol)
    list_reaction_time_theory_nmove.append(T_col_avr(radius_disk=radius_from_length, rho1=0.6, rho2=0.6,
                       diffusion2d=0.1))
areas = [i**2 for i in all_length]
l = list(df1["Lengths"])
areas_table = [i ** 2 for i in l]
df1["Area"] = areas_table


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)
import matplotlib.pyplot as plt
import statistics
from math import sqrt
def plot_confidence_interval(x, values, z=1.96, color='#1f77b4', horizontal_line_width=10):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / sqrt(len(values))

    left = x - horizontal_line_width / 2
    top = mean - confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean + confidence_interval
    plt.plot([x, x], [top, bottom], color=color, linewidth=4)
    plt.plot([left, right], [top, top], color=color, linewidth=4)
    plt.plot([left, right], [bottom, bottom], color=color,  linewidth=4)
    plt.plot(x, mean, 'o', color=color, markersize = 10)

    return mean, confidence_interval


plt.xticks(areas, areas)
plt.plot(areas, list_reaction_time_theory, linewidth = 7, color='#984ea3')
for an_area in areas:
    dd = df1.loc[df1["Area"] == an_area]
    plot_confidence_interval(an_area, dd["Reaction Time"], color='#1f77b4')
plt.hlines(y=28.75, xmin=9, xmax=areas[-1], linewidth=5, color="#d62728")
plt.xlabel(r"Area, $L^2$", size=37)
plt.ylabel(r"Time", size=37)
plt.xticks(size = 37, rotation=45)
plt.yticks(size = 37)
import matplotlib.patches as mpatches
patch1 = mpatches.Patch(color='#d62728', label='Label3')
patch2 = mpatches.Patch(color='#984ea3', label='Label1')
patch3 = mpatches.Patch(color='#1f77b4', label='Theory')
plt.legend(handles=[patch1, patch2, patch3], labels=[r'$T_{C(A+B)\to C(AB)}$', r'$T_{A+B\to AB}$', 'Theory'], prop={'size':35})
plt.locator_params(axis='x', nbins=8)
plt.locator_params(axis='y', nbins=5)
plt.show()