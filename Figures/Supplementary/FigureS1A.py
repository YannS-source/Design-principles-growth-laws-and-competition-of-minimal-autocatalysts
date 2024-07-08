import numpy as np
import ast

################

"""
Create the df from the txt files of the run, and plot
"""

path_serv_move = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS1/Spontaneous_Reaction"
all_run = np.arange(0, 221, 1)
energetic_barrier = 0
all_length = np.arange(4, 17, 1)
diameter_particle = 1

##################
##################
list_all_run = []
list_all_length = []
list_all_times = []

for a_run in all_run:
    for a_length in all_length:
        try:
            file = path_serv_move + '/length={}_diameterMain=1_r_cut=1.1_run={}.txt'.format(a_length,a_run)
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

df1_avr = df1.groupby(["Lengths"])["Reaction Time"].mean()
print("The spontaneous reaction with different L, move:", df1_avr)


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



list_reaction_time_approx = []
for i in areas:
    list_reaction_time_approx.append(i)

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


plt.xticks(areas, areas)
plt.plot(areas, list_reaction_time_theory, linewidth = 5, color='#FFA500')
for an_area in areas:
    dd = df1.loc[df1["Area"] == an_area]
    plot_confidence_interval(an_area, dd["Reaction Time"], horizontal_line_width=0.5)
plt.xlabel(r"Area, $L^2$", size=37)
plt.ylabel(r"$T_{A+B \to AB}$, $\varepsilon_{AB}^+ =0$", size=37)
plt.xticks(size = 35, rotation=45)
plt.yticks(size = 35)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='#1f77b4', label='Label1')
black_patch = mpatches.Patch(color='#FFA500', label='Label2')
plt.legend(handles=[red_patch, black_patch], labels=['Simulations', r'Theory: absorbing trap'], prop={'size':30})
plt.locator_params(axis='x', nbins=7)


