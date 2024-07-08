import numpy as np

# Define the pairwise Wang-Frenkel potential according to equations (S1) & (S2)
def WF_potential_shift(shift=0, eps=1, sigma=1, nu=1, mu=1, rc=1.2, r=1):
    alpha = 2 * nu * (rc / sigma) ** (2 * mu) * ((1 + 2 * nu) / ((2 * nu) * ((rc / sigma) ** (2 * mu) - 1))) ** (
            2 * nu + 1)
    return eps * alpha * ((sigma / (r - shift)) ** (2 * mu) - 1) * ((rc / (r - shift)) ** (2 * mu) - 1) ** (2 * nu)

# Set initial parameters for the potential between particles of different types (e.g., A B), with repulsive forces
the_sigma = 1
repulsive_barrier_epsilon = 10
strong_epsilon = 40
diameter_particle = 1
the_rc = 1.1

# Wang-Frenkel potential
the_WF_potential_start_1_attract = np.array(
    [WF_potential_shift(shift=0, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc,
                        r=i) + repulsive_barrier_epsilon
     for i in np.arange(1, 1 + (the_rc - 1), 1e-3)])


# Inverted Wang-Frenkel potential for the barrier
rmin = the_rc * (3 / (1 + 2 * (the_rc / the_sigma) ** 2)) ** (1 / 2)
the_WF_potential_start_1_repulse = np.array(
    [-WF_potential_shift(shift=the_rc - 1 - (rmin - 1), eps=repulsive_barrier_epsilon, sigma=the_sigma,
                         nu=1, mu=1, rc=the_rc, r=i)
     for i in
     np.arange(1 + (the_rc - 1), 1 + (the_rc - 1) + the_rc - 1 - (rmin - 1),
               1e-3)])

# For graphic reasons, add the zero potential from 1.175 \sigma to 2.
the_WF_0_to_2 = np.array([0 for i in np.arange(1 + (the_rc - 1) + the_rc - 1 - (rmin - 1), 2, 1e-3)])
PM_strong_potential = np.hstack((the_WF_potential_start_1_attract, the_WF_potential_start_1_repulse, the_WF_0_to_2))



# Similarly, set initial parameters for the potential between particles of the same type (AA), without repulsive forces
the_repulsive_epsilon = 0
repulsive_barrier_epsilon = 0
strong_epsilon = 15
the_WF_potential_start_1_attract = np.array(
    [WF_potential_shift(shift=0, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc,
                        r=i) + repulsive_barrier_epsilon
     for i in np.arange(1, 1 + (the_rc - 1), 1e-3)])
the_WF_potential_start_1_repulse = np.array(
    [-WF_potential_shift(shift=the_rc - 1 - (rmin - 1), eps=repulsive_barrier_epsilon, sigma=the_sigma,
                         nu=1, mu=1, rc=the_rc, r=i)
     for i in
     np.arange(1 + (the_rc - 1), 1 + (the_rc - 1) + the_rc - 1 - (rmin - 1), 1e-3)])
PM_strong_potential_waa = np.hstack((the_WF_potential_start_1_attract, the_WF_potential_start_1_repulse, the_WF_0_to_2))


# Set the x axis
x_axis = np.hstack((np.arange(1, 1 + (the_rc - 1), 1e-3), np.arange(1 + (the_rc - 1), 1 + (the_rc - 1) + the_rc - 1 - (rmin - 1),
                1e-3), np.arange(1 + (the_rc - 1) + the_rc - 1 - (rmin - 1),2, 1e-3)))


# Plot both potentials
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(14.5, 10.5)
plt.plot(x_axis, PM_strong_potential, color='#1f77b4',  linewidth=7, label=r"$A+B$ $\rightleftharpoons$ $AB$")
plt.plot(x_axis, PM_strong_potential_waa, color='red',  linewidth=7, label=r"$A+A$ $\rightleftharpoons$ $AA$")
plt.xticks(size = 35)
plt.yticks(size = 35)
desired_xticks = [1.0, 1.03, 1.1, 1.175]
plt.xticks(desired_xticks, [str(x) for x in desired_xticks], size=35)
plt.legend(loc="lower right", prop={'size': 35})
plt.xlim(xmin=1.0, xmax=1.2)
plt.locator_params(axis='y', nbins=5)
plt.show()