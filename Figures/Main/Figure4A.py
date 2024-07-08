import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def estar(a, b, c, x):
    addinh = max(0, 2 * x - b)
    return max(a -np.log(2) + addinh, a, c, 2*x,
               -x+ 2*a -np.log(2) + addinh,
               -2*x + c + 2*a + addinh,
               -x+c+a + np.log(2))
def tcycle(a, b, c, x):
    addinh = 2 * x - b
    G1, G2, G3, G4 = a -np.log(2), a, c, 2*x
    G5, G6, G7 = -x + 2*a-np.log(2), -2*x + c + 2*a, -x+c+a + np.log(2)
    G8, G9, G10 = a -np.log(2) + addinh, -x+ 2*a -np.log(2) + addinh, -2*x + c + 2*a + addinh
    # print(G1, G2, G3, G4)
    # print(G5, G6, G7)
    # print(G8, G9, G10)
    return np.exp(G1) + np.exp(G2) + np.exp(G3) + np.exp(G4) + np.exp(G5) +\
           np.exp(G6) + np.exp(G7) + np.exp(G8) + np.exp(G9) + np.exp(G10)
def regime(a, c, x, b=np.inf, spo=False):
    addinh = max(0, 2 * x - b)
    es = estar(a, b, c, x)
    #if spo and es > 2*a + c:
    if spo and np.exp(2*a + c) < tcycle(a, b, c, x):
        return 10
    if es == a -np.log(2) + addinh:
        if addinh == 0:
            return 0
        else:
            return 7
    if es == a:
        return 1
    if es == c:
        return 2
    if es == 2 * x:
        return 3
    if es == -x+ 2*a -np.log(2) + addinh:
        if addinh == 0:
            return 4
        else:
            return 8
    if es == -2*x + c + 2*a + addinh:
        if addinh == 0:
            return 5
        else:
            return 9
    if es == -x+c+a + np.log(2):
        return 6
def xopt(a, c, xrange=np.arange(0, 1.55, .01), b=np.inf):
    xo = min(xrange)
    eo = estar(a, b, c, xo)
    for x in xrange:
        if estar(a, b, c, x) < eo:
            xo, eo = x, estar(a, b, c, x)
    return xo
def xopt_total(a, c, xrange=np.arange(0, 1.55, .01), b=np.inf):
    xo = min(xrange)
    eo = tcycle(a, b, c, xo)
    for x in xrange:
        if estar(a, b, c, x) < eo:
            xo, eo = x, estar(a, b, c, x)
    return xo

"""Code"""

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


# Plot setup
fig, ax = plt.subplots(figsize=(14, 10))

le_c = 10
window= 1
le_b = 3/window # Concentration [AB]
arange = np.arange(0, 10/window, 0.1/window)
xrange = np.arange(0, 10/window, 0.1/window)
cols = [regime_labels[i][0] for i in range(len(regime_labels))]

used_regimes = set()  # To track which regimes are used

for a in arange:
    for x in xrange:
        reg = regime(a, le_c, x, b=le_b, spo=True)
        ax.scatter(a, x, c=cols[reg], s=100)
        used_regimes.add(reg)

used_regimes.add(11)
used_regimes = sorted(used_regimes)
# For aestethics want the spontaneous reaction to come first
try:
    used_regimes.pop(used_regimes.index(10))
    used_regimes.insert(0, 10)
except ValueError:
    pass
xopt_val = [xopt(a, le_c, xrange=xrange, b=le_b) for a in arange]
xopt_tot = [xopt_total(a, le_c, xrange=xrange, b=le_b) for a in arange]

# a small problem of precision that makes the figure ugly, remove it simply manually:
new_xopt_tot =[]
for i in range(len(xopt_tot)):
    if i < 30:
        new_xopt_tot.append(xopt_tot[i])
    else:
        new_xopt_tot.append(1.5)



"""Plot"""

ax.plot(arange, new_xopt_tot, color='red', lw=8)

# Generate legend patches
patches = [mpatches.Patch(color=regime_labels[reg][0], label=regime_labels[reg][1]) for reg in used_regimes[:-1]]
patches.append(mpatches.Patch(color=regime_labels[11][0], label=""))
plt.legend(handles=patches, loc='upper left', prop={'size': 28}, framealpha=0.8)

plt.xticks(size=39)
plt.yticks(size=39)
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.xlim(xmin = 0.0, xmax = 10.0)
plt.ylim(ymin = 0.0, ymax = 10.0)

for spine in ax.spines.values():
    spine.set_visible(False)
plt.tight_layout()
plt.show()