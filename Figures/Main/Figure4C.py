import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def cycling_time_w_product_inhibition(xA, xC, k1, k_1, k3, k4, k_4):
    pIp, pIm = k4, k_4 * xC
    p1p, p1m = 2 * xA * k1, k_1
    p2p, p2m = xA * k1, 2 * k_1
    p3p = k3
    p4p = k4
    return 1 / p1p + 1 / p2p + 1 / p3p + 1 / p4p + p1m / (p1p * p2p) + (p1m * p2m) / (p1p * p2p * p3p) + p2m / (
            p2p * p3p) + (pIm / pIp) * (1 / p1p + p1m / (p1p * p2p) + (p1m * p2m) / (p1p * p2p * p3p))

file = "/Users/yannsakref/Dev/Dev_PhD/Paper3/Figures/Main/Fig4/start_0_AA_omega_c_10_delta_d_0_ABtot_3.csv"
# The Markov model detailed in the method section is numerically integrated for different paramaters until [AB]_{tot} = 3
# At this point, the rate of the spontaneous reaction k_{sp}, the reaction order n, the reaction rate k,
# the limiting barrier, and the concentration of the different species are reported in the previous file.
# The following code plot the figure 4C, the reaction order n, based on this file.

df = pd.read_csv(file)
le_b = 3

maximal_xc = 10.1
all_wc = [1]
all_detlad = [0]


"""DO n First"""
df = df[(df["x/c"] >= 0) & (df["x/c"] < maximal_xc)]
df = df[(df["a/c"] >= 0) & (df["a/c"] < maximal_xc)]
# Condition where k_sp is greater than k
condition = df['k_sp'] > df['k']
# Set n to 0 where the condition is True
df.loc[condition, 'n'] = 0
df.loc[df['n'] == 0, 'n'] = -1
df.loc[condition, 'limiting barrier'] = "b11"


list_n = df["n"]
list_n = [0.999 if value >= 1 else value for value in list_n]
df["n"] = list_n
df.loc[df['n'] > 1, 'n'] = 1
dfs = df.pivot(index="x/c", columns='a/c', values='n')
dfs = dfs.fillna(0.00)

"""Do k Now"""
dfk = pd.read_csv(file)
dfk = dfk[(dfk["x/c"] >= 0) & (dfk["x/c"] < maximal_xc)]
dfk = dfk[(dfk["a/c"] >= 0) & (dfk["a/c"] < maximal_xc)]


new_list_a_c = []
list_all_weak_epsilon = []
list_all_relative_k = []
column_to_plot = 'k' # could be "k", 'rate_Tc', "k_AB01", "rate_Tc_AB01"

list_concentrations_a_c = df["a/c"].unique()
list_concentrations_a_c = [i*all_wc[0] for i in list_concentrations_a_c]
all_weak_epsilon_x_c = df["x/c"].unique()
all_weak_epsilon_x_c = [i*all_wc[0] for i in all_weak_epsilon_x_c]

for a_c in list_concentrations_a_c:
    print(a_c)
    dfx = dfk[(dfk["a/c"] == a_c)]
    print(dfx)
    dfx[column_to_plot] = dfx[column_to_plot] / np.nanmax(dfx[column_to_plot])
    for a_we in range(len(all_weak_epsilon_x_c)):
        print(all_weak_epsilon_x_c[a_we])
        new_list_a_c.append(a_c)
        list_all_weak_epsilon.append(all_weak_epsilon_x_c[a_we])
        if pd.isna(list(dfx[column_to_plot])[a_we]):
            list_all_relative_k.append(0)
        else:
            list_all_relative_k.append(list(dfx[column_to_plot])[a_we])


dfk = pd.DataFrame(list(zip(list_all_weak_epsilon, new_list_a_c, list_all_relative_k)),
               columns =["x/c", '$a/c$', column_to_plot])

dfsk = dfk.pivot(index="x/c", columns='$a/c$', values=column_to_plot)






###########
"""Plot"""

import seaborn as sns
from scipy import ndimage


fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)

smooth_scale = 1
z = ndimage.zoom(dfs.to_numpy(), smooth_scale)
z = np.clip(z, a_min=None, a_max=1)  # ensure interpolated values do not exceed 1

# Create a zoomed and clipped array for dfsk
zk = ndimage.zoom(dfsk.to_numpy(), smooth_scale)
zk = np.clip(zk, a_min=None, a_max=1)  # ensure interpolated values do not exceed 1


# When plotting, use a colormap with 'dimgrey' for the negative value
cmap = plt.cm.Blues  # Or any other colormap
# cmap = plt.cm.Reds  # Or any other colormap
cmap.set_under("#A6A6A6")  # Set color for under-range values
# cntr = ax.contourf(np.linspace(0, len(dfs.columns), len(dfs.columns) * smooth_scale),
#                    np.linspace(0, len(dfs.index), len(dfs.index) * smooth_scale),
#                    z, levels=np.arange(0.5, 1.099, 0.1), cmap='Oranges_r')

cntr = ax.contourf(np.linspace(0, len(dfs.columns), len(dfs.columns) * smooth_scale),
                   np.linspace(0, len(dfs.index), len(dfs.index) * smooth_scale),
                   z, levels=np.arange(0.5, 1.099, 0.1), cmap=cmap, extend='min')

# Adding contour lines for "n" = 0.5
ax.contour(np.linspace(0, len(dfs.columns), len(dfs.columns) * smooth_scale),
           np.linspace(0, len(dfs.index), len(dfs.index) * smooth_scale),
           z, levels=[0.9], linewidths=7, colors = "black")#colors='chartreuse')

# cntr = ax.contourf(np.linspace(0, len(dfs.columns), len(dfs.columns) * smooth_scale),
#                    np.linspace(0, len(dfs.index), len(dfs.index) * smooth_scale),
#                    z, levels=[0, 0.1], cmap='Gray')
# Adding contour lines for "n" = 0.5 for dfsk
# ax.contour(np.linspace(0, len(dfsk.columns), len(dfsk.columns) * smooth_scale),
#            np.linspace(0, len(dfsk.index), len(dfsk.index) * smooth_scale),
#            zk, levels=[0.99], colors='red', linewidths=7)  # Set linewidths to 3 for thickness

# Find the maximum value in each 'OmegaC' column and its corresponding 'x/c' index
max_vals = dfsk.idxmax()  # This gives the index (x/c value) of the maximum in each column
max_ac = dfsk.columns  # OmegaC values
max_xc_values = [max_vals[omega] for omega in max_ac]  # Extracting the 'x/c' values corresponding to each max

# Convert OmegaC values and corresponding x/c indices to plottable lists
omega_c_list = max_ac.astype(float).tolist()  # Ensure OmegaC values are float for plotting
xc_list = max_xc_values  # These are the 'x/c' indices (as floats, if not, convert them)
# Correction: Ensure that the x and y coordinates for the line plot are mapped correctly
# Since your contour plot uses np.linspace() to define its axes, align your line plot accordingly
# Plot the line connecting the maximum points
# Convert xc_list (y-axis) into indices to match the linspace used in the contour plot
xc_indices = [np.where(dfs.index == xc_val)[0][0] for xc_val in xc_list]  # Convert x/c values to their indices in dfs.index
omega_indices = [np.where(dfs.columns == omega_val)[0][0] for omega_val in omega_c_list]  # Similar conversion for OmegaC, if necessary

# ax.plot(omega_indices, xc_indices, color='red', linewidth=5, marker='o', linestyle='-', markersize=8)
ax.plot(omega_indices, xc_indices, color='red', linewidth=8, linestyle='-', alpha=1)

ax = sns.heatmap(dfs, annot=False, alpha=0, cbar=False, ax=ax)
cbar = plt.colorbar(cntr, ax=ax)
tick_font_size = 35
cbar.ax.tick_params(labelsize=tick_font_size)

import matplotlib.patches as mpatches
patches = []
patches.append(mpatches.Patch(color="#A6A6A6", label=r'$T_{\rm cycle} < T_{sp}$'))
patches.append(mpatches.Patch(color='red', label=r''))
plt.legend(handles=patches, loc='upper left', prop={'size': 27}, framealpha=0.9)

# Set x-axis ticks and labels
xticks = np.arange(len(dfs.columns))
xticklabels = dfs.columns.values
# Filter to keep only x-ticks at every 0.5 interval
filtered_xticks = xticks[xticklabels % 2 == 0]
filtered_xticklabels = xticklabels[xticklabels % 2 == 0]
filtered_xticklabels = filtered_xticklabels.astype(int)
# Setting the filtered ticks and labels to the axis
ax.set_xticks(filtered_xticks)
ax.set_xticklabels(filtered_xticklabels)
# Set y-axis ticks and labels
yticks = np.arange(len(dfs.index))
yticklabels = dfs.index.values
print(yticklabels)
# Filter to keep only y-ticks at every 0.5 interval
filtered_yticks = yticks[np.isclose(yticklabels % 1, 0)]
filtered_yticklabels = yticklabels[np.isclose(yticklabels % 1, 0)]
print(filtered_yticklabels)
filtered_yticklabels = filtered_yticklabels.astype(int)
# Setting the filtered ticks and labels to the axis
ax.set_yticks(filtered_yticks)
ax.set_yticklabels(filtered_yticklabels)
ax.invert_yaxis()



plt.locator_params(axis='x', nbins=10)
plt.locator_params(axis='y', nbins=10)
plt.xticks(size = 35)
plt.yticks(size = 35)
plt.xticks(rotation=360)

ax = plt.gca()
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

components = file.split("/")
target_component = components[-1]
final_component = target_component.replace(".csv", "")

plt.gca().set_xlabel('')
plt.gca().set_ylabel('')
# plt.savefig(
#     "/Users/yannsakref/Dev/Dev_PhD/Paper3/Figures/Main/Fig4/wc_10_real_n_b={}.pdf".format(le_b),
#     format="pdf", bbox_inches="tight")
# plt.close()
plt.show()