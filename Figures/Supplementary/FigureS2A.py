import numpy as np
import pandas as pd

df = pd.read_csv("/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/GitHub/Data/FigureS2/Release_A_AB_to_13_e5.csv")

df_onlyA = df.loc[df["Type"] == 'Only_A_']
df_AB = df.loc[df["Type"] == 'Release_AB_']

df_mean_Tr_over_epsilon = df_onlyA.groupby(['Weak Epsilon'])["Release time"].mean()
df_mean_Tr2_over_epsilon = df_AB.groupby(['Weak Epsilon'])["Release time"].mean()

print(df_mean_Tr_over_epsilon)
print(df_mean_Tr2_over_epsilon)

all_epsilon = np.arange(0, 14, 1)

v = [0.823842 * np.exp((i-6)/1.3) for i in all_epsilon]


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


plt.xticks(all_epsilon, all_epsilon)
plt.plot(all_epsilon, v, linewidth = 5, color='#FFA500', label=r'$T_{A\cdot AB +B\to A+B+AB}(\epsilon_{AA} = 6)e^{(\epsilon_{AA} - 6)}$')
for an_eps in all_epsilon:
    dd = df_onlyA.loc[df_onlyA["Weak Epsilon"] == an_eps]
    plot_confidence_interval(an_eps, dd["Release time"] )
plt.xlabel(r"$\epsilon_{AA}^+$", size=37)
plt.ylabel(r'$T_{A\cdot AB +B\to A+B+AB}$', size=37)
plt.yscale("log")
plt.xticks(size = 35)
plt.yticks(size = 35)
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='#1f77b4', label='Label1')
black_patch = mpatches.Patch(color='#FFA500', label='Label2')
plt.legend(handles=[red_patch, black_patch], labels=['Simulations', r'$T_{A\cdot AB\to A+AB}(\epsilon_{AA} = 6)e^{(\epsilon_{AA} - 6)}$'], prop={'size':30})
plt.locator_params(axis='x', nbins=7)
plt.show()






