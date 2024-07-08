import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from scipy import ndimage

file = "/Users/yannsakref/Dropbox/Devs/Dev/Dev_PhD/Paper3/Data/Figure4/start_0_AA_omega_c_1_delta_d_0_ABtot_3.csv"
# The Markov model detailed in the method section is numerically integrated for different paramaters until [AB]_{tot} = 3
# At this point, the rate of the spontaneous reaction k_{sp}, the reaction order n, the reaction rate k,
# the limiting barrier, and the concentration of the different species are reported in the previous file.
# The following code plot the figure 4B, the limiting barrier, based on this file.

df = pd.read_csv(file)
all_wc = [10]
maximal_xc = 10.1
all_weak_epsilon_x_c = df["x/c"].unique()
all_weak_epsilon_x_c = [i*1 for i in all_weak_epsilon_x_c]
all_weak_epsilon_x_c = [x for x in all_weak_epsilon_x_c if x < maximal_xc]
le_b = 3

"""First, compute the limiting barriers"""
# Convert the 'limiting barrier' column to numeric values
df = df[(df["x/c"] >= 0) & (df["x/c"] < maximal_xc)]
df = df[(df["a/c"] >= 0) & (df["a/c"] < maximal_xc)]
condition = df['k_sp'] > df['k']
df.loc[condition, 'limiting barrier'] = "b11"
df['limiting barrier'] = df['limiting barrier'].apply(lambda x: int(x[1:])) -1
# Pivot the DataFrame
pivot_df = df.pivot_table(index='x/c', columns='a/c', values='limiting barrier', aggfunc=np.mean)
print(pivot_df)


"""Do k Now"""
dfk = pd.read_csv(file)
dfk = dfk[(dfk["x/c"] >= 0) & (dfk["x/c"] < maximal_xc)]
dfk = dfk[(dfk["a/c"] >= 0) & (dfk["a/c"] < maximal_xc)]


new_list_a_c = []
list_all_weak_epsilon = []
list_all_relative_k = []
column_to_plot = 'k' # could be "k", 'rate_Tc', "k_AB01", "rate_Tc_AB01"

list_concentrations_a_c = df["a/c"].unique()
list_concentrations_a_c = [i*1 for i in list_concentrations_a_c]
all_weak_epsilon_x_c = df["x/c"].unique()
all_weak_epsilon_x_c = [i*1 for i in all_weak_epsilon_x_c]

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


"""Plotting"""

# Color and label mapping for each regime

regime_labels = {
    0: ("w", u'$G_1^{\u2021} - G_1$'),
    1: ("green", u'$G_2^{\u2021} - G_2$'),
    2: ("skyblue", u'$G_3^{\u2021} - G_3$'),
    3: ("#FFA500", u'$G_4^{\u2021} - G_4$'),
    4: ("#7EFF01", u'$G_2^{\u2021} - G_1$'),
    5: ("lightgreen", u'$G_3^{\u2021} - G_1$'),
    6: ("#208FFF", u'$G_3^{\u2021} - G_2$'),
    7: ("firebrick", u'$G_1^{\u2021} - G_0$'),
    8: ("green", u'$G_2^{\u2021} - G_0$'),
    9: ("teal", u'$G_3^{\u2021} - G_0$'),
    10: ("#A6A6A6", r'$T_{\rm cycle} < T_{sp}$'),
11: ("red", r'$\hat \epsilon_{AA}^-$')
}

used_regimes = set()  # To track which regimes are used

# Custom colors and labels
cols = [list(regime_labels.values())[i][0] for i in range(len(list(regime_labels.values())))]


# Convert colors to RGBA
cols_rgba = [mcolors.to_rgba(color) for color in cols]

# Map the values in pivot_df_filled to RGBA colors
color_matrix_rgba = np.array([[cols_rgba[val] for val in row] for row in pivot_df.values])
# Track the used regimes
used_regimes = set(pivot_df.values.ravel())
used_regimes.add(11)
used_regimes = sorted(used_regimes)

try:
    used_regimes.pop(used_regimes.index(10))
    used_regimes.insert(0, 10)
except ValueError:
    pass


"""Plot"""

# Plot the heatmap
fig, ax = plt.subplots(figsize=(14, 10))
cax = ax.imshow(color_matrix_rgba, aspect='auto')
xticks = np.arange(len(pivot_df.columns))
xticklabels = pivot_df.columns.values
# Filter to keep only x-ticks at every 0.5 interval
filtered_xticks = xticks[xticklabels % 2 == 0]
filtered_xticklabels = xticklabels[xticklabels % 2 == 0]
filtered_xticklabels = filtered_xticklabels.astype(int)
# Setting the filtered ticks and labels to the axis
ax.set_xticks(filtered_xticks)
ax.set_xticklabels(filtered_xticklabels)


# Set y-axis ticks and labels
yticks = np.arange(len(pivot_df.index))
yticklabels = pivot_df.index.values
print(yticklabels)
# Filter to keep only y-ticks at every 0.5 interval
filtered_yticks = yticks[np.isclose(yticklabels % 1, 0)]
filtered_yticklabels = yticklabels[np.isclose(yticklabels % 1, 0)]
print(filtered_yticklabels)
filtered_yticklabels = filtered_yticklabels.astype(int)
# Setting the filtered ticks and labels to the axis
ax.set_yticks(filtered_yticks)
ax.set_yticklabels(filtered_yticklabels)

# Interpolate dfsk for smoother contours
smooth_scale = 10  # Change this value as needed for smoother contours
zk = ndimage.zoom(dfsk.to_numpy(), smooth_scale)

# Adjust the linspace to match the heatmap's axes
x_contour = np.linspace(0, len(pivot_df.columns) - 1, len(dfsk.columns) * smooth_scale)
y_contour = np.linspace(0, len(pivot_df.index) - 1, len(dfsk.index) * smooth_scale)


# Find the maximum value in each 'OmegaC' column and its corresponding 'x/c' index
max_vals = dfsk.idxmax()  # This gives the index (x/c value) of the maximum in each column
max_omegas = dfsk.columns  # OmegaC values
max_xc_values = [max_vals[omega] for omega in max_omegas]  # Extracting the 'x/c' values corresponding to each max

# Convert OmegaC values and corresponding x/c indices to plottable lists
omega_c_list = max_omegas.astype(float).tolist()  # Ensure OmegaC values are float for plotting
xc_list = max_xc_values  # These are the 'x/c' indices (as floats, if not, convert them)
# Correction: Ensure that the x and y coordinates for the line plot are mapped correctly
# Since the contour plot uses np.linspace() to define its axes, align your line plot accordingly
# Plot the line connecting the maximum points
# Convert xc_list (y-axis) into indices to match the linspace used in the contour plot
xc_indices = [np.where(pivot_df.index == xc_val)[0][0] for xc_val in xc_list]  # Convert x/c values to their indices in dfs.index
omega_indices = [np.where(pivot_df.columns == omega_val)[0][0] for omega_val in omega_c_list]  # Similar conversion for OmegaC, if necessary

ax.plot(omega_indices, xc_indices, color='red', linewidth=8, linestyle='-', alpha=1)
# Generate legend patches
patches = [mpatches.Patch(color=regime_labels[reg][0], label=regime_labels[reg][1]) for reg in used_regimes[:-1]]
patches.append(mpatches.Patch(color=regime_labels[11][0], label=""))


for spine in ax.spines.values():
    spine.set_visible(True)


plt.legend(handles=patches, loc='upper left', prop={'size': 28}, framealpha=0.8)

plt.xticks(size=35)
plt.yticks(size=35)
plt.locator_params(axis='x', nbins=10)
plt.locator_params(axis='y', nbins=10)
ax.invert_yaxis()
plt.show()
