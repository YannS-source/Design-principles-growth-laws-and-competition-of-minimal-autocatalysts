import pandas as pd
import numpy as np


df = pd.read_csv("/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/Data/Figure5/wc_10_A04.5399929762484854e-05_wAA8_wBB8_wCC8_wDD14_wEE20.csv")
# df with the steady-state concentrations of the different species by themselves or when competiting
# for a comon resource

all_logphis = np.arange(12, 28, 0.1)

concentration_AB_alone = df["loneAB"]
concentration_AC_alone = df["loneAC"]
concentration_AD_alone = df["loneAD"]

concentrations_AB = df["AB"]
concentrations_AC = df["AC"]
concentrations_AD = df["AD"]

all_logphis = [np.exp(i-10) for i in all_logphis]
log_all_logphis = [np.log(i) for i in all_logphis]
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)


plt.plot(log_all_logphis[:-1], concentrations_AC[:-1],  linewidth=7, color ="blue", label=r"$AB$, $\epsilon_{AB}^+=10$, $\epsilon_{BB}^-=8$")
plt.plot(log_all_logphis[:-1], concentrations_AD[:-1],  linewidth=7, color ="#8c564b", label="$AC$, $\epsilon_{AC}^+=10$, $\epsilon_{CC}^-=12$")
plt.plot(log_all_logphis[:-1], concentrations_AB[:-1],  linewidth=7, color ="green", label="$AD$, $\epsilon_{AD}^+=12$, $\epsilon_{DD}^-=8$")
plt.plot(log_all_logphis[:-1], concentration_AC_alone[:-1],  linewidth=7, linestyle="dotted", color ="blue")
plt.plot(log_all_logphis[:-1], concentration_AD_alone[:-1],  linewidth=7, linestyle="dotted", color ="#8c564b")
plt.plot(log_all_logphis[:-1], concentration_AB_alone[:-1],  linewidth=7, linestyle="dotted", color ="green")
plt.locator_params(axis='x', nbins=4)

plt.xlabel(r'$\tau/T_{C(A+B)\to CAB}$', size=40)
plt.ylabel(r'Steady-State Concentration/$[A_0]$', size=40)
plt.xticks(size = 37)
plt.yticks(size = 37)
plt.legend(prop={'size':35})

# Adjusting layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.85, 1])


plt.show()