import numpy as np
import pandas as pd

### THIS IS THE SIMULATIONS
df = pd.read_csv("/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS2/Release_A_AB_to_13_e5.csv")

df_onlyA = df.loc[df["Type"] == 'Only_A_']
df_AB = df.loc[df["Type"] == 'Release_AB_']
df_mean_Tr_over_epsilon = df_onlyA.groupby(['Weak Epsilon'])["Release time"].mean()
df_mean_Tr2_over_epsilon = df_AB.groupby(['Weak Epsilon'])["Release time"].mean()
we = list(df_mean_Tr_over_epsilon.index)

## THIS IS THE THEORY
the_sum = np.linspace(0, 10000, 10001)
def min_t_2_t3(radius, D, alpha, Tr):
    h = 1/Tr
    f1 = (D) / ((2 * np.pi * radius - alpha * 2 * np.pi * radius) ** 2)
    temps_min_fp = 0
    for n in the_sum:
        prefactor = 8
        other = ((np.pi * (2 * n + 1)) ** 2) * (f1 * ((np.pi * (2 * n + 1)) ** 2) + 2 * h)
        temps_min_fp += prefactor / other
    return temps_min_fp
def t_2_t3(radius, D, alpha, Tr):
    h = 1 / Tr
    f1 = (D) / ((2 * np.pi * radius - alpha * 2 * np.pi * radius) ** 2)
    proba_retour_aller_fp = 0
    for n in the_sum:
        prefactor = 16 * h
        other = ((np.pi * (2 * n + 1)) ** 2) * (f1 * ((np.pi * (2 * n + 1)) ** 2) + 2 * h)
        proba_retour_aller_fp += prefactor / other
    return proba_retour_aller_fp
def min_t_4_t5(radius, D, Tr, Gamma):
    h = 1 / Tr
    f2 = (D) / ((np.pi * radius - Gamma) ** 2)
    temps_min_sp = 0
    for n in the_sum:
        prefactor = 4 * np.sin((x0 * np.pi * (2 * n + 1)) / (np.pi*radius - Gamma))
        other = (np.pi * (2 * n + 1)) * (f2 * ((np.pi * (2 * n + 1)) ** 2) + h)
        temps_min_sp += prefactor / other
    return temps_min_sp
def t_4_t5(radius, D, Tr, Gamma):
    h = 1 / Tr
    f2 = (D) / ((np.pi * radius - Gamma) ** 2)
    proba_retour_aller_sp = 0
    for n in the_sum:
        prefactor = 4 * (np.pi * (2 * n + 1)) * f2 * np.sin((x0 * np.pi * (2 * n + 1)) / (np.pi*radius - Gamma))
        other = (f2 * ((np.pi * (2 * n + 1)) ** 2) + h)
        proba_retour_aller_sp += prefactor / other
    return proba_retour_aller_sp

x0 = 0.093
a_gamma = 0.1
a_radius = 0.5
a_bulk_diffusion = 0.1
a_rota_diffusion = 0.1
all_Tr = list(df_mean_Tr_over_epsilon)
all_beta_Tcol = [10]


time_exp = []
time_normal = []
time_spontaneous_reaction = []
list_of_big_radii = []
list_of_Tr = []
list_of_radius = []
list_of_rota_diffusion = []
list_of_Tsp_Tc = []
list_of_beta_tcol = []
proba_bind_over_release = []
proba_bind_over_release1 = []
first_part = []
second_part = []
ratio_both_part = []
for a_Tr in all_Tr:
    h = 1/ a_Tr

    temps_min_sp = min_t_4_t5(radius=a_radius, D=a_rota_diffusion, Tr=a_Tr, Gamma=a_gamma)
    proba_retour_aller_sp = t_4_t5(radius=a_radius, D=a_rota_diffusion, Tr=a_Tr, Gamma=a_gamma)


    A2 = np.array([[1, -1],
                   [-proba_retour_aller_sp, 1]])

    b2 = np.array([1 / (2 * h),
                   temps_min_sp])

    x2 = np.linalg.solve(A2, b2)

    list_of_Tr.append(a_Tr)
    list_of_rota_diffusion.append(a_rota_diffusion)
    second_part.append(x2[0])

list_epsilon = list(df_AB["Weak Epsilon"])
dico_list_epsilon = {}
for i in np.arange(0, 14, 1):
    dico_list_epsilon[i] = list_of_Tr[i]
release_A_in_AB = []
for i in list_epsilon:
    if i in np.arange(0, 14, 1):
        release_A_in_AB.append(dico_list_epsilon[i])
df_AB["Release_A"] = release_A_in_AB
v = [1 * np.exp(2*((i-6)/1.3)) for i in we]



## Plotting

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

plt.xticks(we, we)
plt.plot(we, second_part,  linewidth=5,color='#FFA500', label = "Theory")
plt.plot(we, v, linewidth = 5, color="#2ca02c", label=r'$T_{sp}(h_{AB}^- = 6)e^{h_{AB}^- - 6}$')
print("Second part:", second_part)

for a_we in we:
    dd = df_AB.loc[df_AB["Weak Epsilon"] == a_we]
    plot_confidence_interval(a_we, dd["Release time"], horizontal_line_width=0.1)


plt.xlabel(r'$\epsilon_{AA}$', size=37)
plt.ylabel(r'$T_{AB{:}AB \to AB+AB}$', size=37)
plt.yscale("log")
plt.xticks(size = 35)
plt.yticks(size = 35)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='#1f77b4', label='Label1')
black_patch = mpatches.Patch(color='#FFA500', label='Label2')
black_patch2 = mpatches.Patch(color='#2ca02c', label='Label2')
plt.legend(handles=[red_patch, black_patch, black_patch2], labels=['Simulations', r'Theory: $1D$ Diffusion', r"$e^{2((\epsilon_{AA}-6)/1.3))}$"], prop={'size':30})
plt.locator_params(axis='x', nbins=7)
plt.show()
