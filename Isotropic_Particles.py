# Import necessary libraries
import gsd.hoomd
import hoomd
import numpy as np
import os

print(hoomd.version.version)

# Function to compute position based on angle
def compute_position_based_on_angle(diameter, angle):
    return [(diameter / 2) * np.cos(angle), (diameter / 2) * np.sin(angle), 0]
# Function to check if AB and CD have formed in a snapshot
def check_two_dimers_C1C2_AB_one_step(snapshot_position, list_2dimensions, rcut=1.1):
        ## Retun true if two dimers have formed:
        # [d(A-B) < rcut, d(C-D) < rcut] AND d(A-C) > rcut, d(B-D) > rcut
        posC1 = snapshot_position[0][:-1]
        posC2 = snapshot_position[1][:-1]
        posA = snapshot_position[2][:-1]
        posB = snapshot_position[3][:-1]
        distanceC1C2 = distance(posC1, posC2, list_2dimensions)
        if distanceC1C2 > rcut:
            print("WARNING: The catalyst has dissociated")
        distanceAB = distance(posA, posB, list_2dimensions)
        distanceC1A = distance(posC1, posA, list_2dimensions)
        distanceC2B = distance(posC2, posB, list_2dimensions)
        if distanceC1A > rcut and distanceC2B > rcut and distanceAB < rcut:
            return True
# Function to calculate distance between two points considering periodic boundary conditions
def distance(p1, p2, list_2dimensions):
    total = 0
    for i, (a, b) in enumerate(zip(p1, p2)):
        delta = abs(b - a)
        if delta > list_2dimensions[i] - delta:
            delta = list_2dimensions[i] - delta
        total += delta ** 2
    return total ** 0.5

def Model1(folder_to_store='/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Model1/test',
                                   length_box=5, diameter_particle=1,
                                   list_positions=[[0, 0, 0], [1.05, 0, 0], [0, 1.05, 0], [1.05, 1.05, 0]],
                                   list_orientation=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]],
                                   weak_epsilon=0.5, strong_epsilon=40, repulsive_barrier_epsilon=1, the_dt=1e-5, nb_steps=1e5, length_step = 100,
                                   the_recoding_period=100, run=0, drag_coeff=10, the_kT=1,
                                   stop_at_event=True, saving_trajectory=True, device = "cpu", the_rc = 1.1, the_repulsive_rc = 1.1, number_beyond = 1.1):

    """
    Input:
    length_box: length of the square box, with periodic boundary conditions.
    diameter_particle: diameter of A, B, C, and D.
    list_position: list of the particle's position at initialization, [[A], [B], [C], [D]],
    where [i] means [x_i, y_i, z_i].
    list_orientation: list of the particle's position at start. Of no consequence for isotropic particles.
    weak_epsilon: interaction strength between A-C and B-D
    strong_epsilon: interaction strength between A-B and C-D
    repulsive_barrier_epsilon: interaction barrier between A and B.
    the_dt: time step
    nb_steps: number of steps after which the simulation stops.
    length_step: length of the steps
    the_recording_period: in case saving_trajectory=True, the number of step between saving snapshots.
    run: random seed used.
    drag_coeff: the drag coefficient. Generally 10, such that the diffusion coefficient is 0.1.
    the_kT: temperature, 1 by default.
    stop_at_event: stop at a specific event (here, when AB and CD) have formed and dissociated from each other.
    device: 'cpu' or 'gpu'.
    the_rc: distance between 2 particles beyond which the attractive potentials are not 'felt'.
    the_repulsive_rc: distance between 2 particles beyond which the repulsive potentials are not 'felt'.
    number_beyond: the distance between 2 particles where it is considered they do not interact anymore -- this is useful in the
    case where we want to enforce a specific distance between particles to stop the simulation, e.g., that the two
    dimer AB and CD are far enough from each other (typically further than the_rc).

    Returns 3 arguments if stop_at_event==True, none otherwise:
    broke: 1 if the catalyst CD has broken during the simulation, 0 otherwise.
    state_3: time at which the A and B have first formed a bound.
    time: final time, when the two dimers AB and CD have formed.
    """
    ##############################################
    ###            Initialize snapshot        ###
    #############################################
    L = length_box
    num_bodies = 4 # A, B, C, and D
    snapshot = gsd.hoomd.Snapshot()
    snapshot.particles.N = num_bodies
    snapshot.particles.position = list_positions
    snapshot.particles.diameter = [diameter_particle] * num_bodies
    snapshot.particles.moment_inertia = [[0, 0, 1]] * num_bodies
    snapshot.particles.types = ["C1", "C2", "A", "B"] # C1 = C, C2 = D
    snapshot.configuration.box = [L, L, 0, 0, 0, 0]
    snapshot.particles.orientation = list_orientation
    snapshot.particles.typeid = [0, 1, 2, 3]

    with gsd.hoomd.open(name=folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon), mode='wb') as f:
        f.append(snapshot)

    integrator = hoomd.md.Integrator(dt=the_dt)
    cell = hoomd.md.nlist.Cell(buffer=0)

    ##############################################
    ###            Specify interactions        ###
    #############################################

    # Implement the Wang-Frankel potential and related force, according to Eq. (S1), (S2), and (S3).
    def WF_potential_shift(shift=0, eps=1, sigma=1, nu=1, mu=1, rc=1.2, r=1):
        alpha = 2 * nu * (rc / sigma) ** (2 * mu) * ((1 + 2 * nu) / ((2 * nu) * ((rc / sigma) ** (2 * mu) - 1))) ** (
                2 * nu + 1)
        return eps * alpha * ((sigma / (r - shift)) ** (2 * mu) - 1) * ((rc / (r - shift)) ** (2 * mu) - 1) ** (2 * nu)
    def WF_force_shift(shift=0, eps=1, sigma=1, nu=1, mu=1, rc=1.2, r=1):
        alpha = 2 * nu * (rc / sigma) ** (2 * mu) * ((1 + 2 * nu) / ((2 * nu) * ((rc / sigma) ** (2 * mu) - 1))) ** (
                2 * nu + 1)
        CC = (rc / (r - shift)) ** (2 * mu)
        SS = (sigma / (r - shift)) ** (2 * mu)
        return (2 / (r - shift)) * eps * alpha * mu * (CC - 1) ** (2 * nu - 1) * (2 * nu * (SS - 1) * CC + SS * (CC - 1))
    the_sigma = 1
    the_repulsive_epsilon = 1

    # Weak attractive potential without repulsive forces, for A-C and B-D
    PM_weak_potential = np.array(
        [WF_potential_shift(shift=0, eps=weak_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange(1, the_rc, 1e-2)])
    PM_weak_force = np.array(
        [WF_force_shift(shift=0, eps=weak_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange(1, the_rc, 1e-2)])

    # Repulsive potential for excluded volume (all particles)
    WF_reject_particle_potential = np.array(
        [WF_potential_shift(shift=(diameter_particle-1), eps=the_repulsive_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange((diameter_particle-1) + 0.01, (diameter_particle-1) + 1.0, 1e-2)])
    WF_reject_particle_force = np.array(
        [WF_force_shift(shift=(diameter_particle-1), eps=the_repulsive_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange((diameter_particle-1) + 0.01, (diameter_particle-1) + 1.0, 1e-2)])

    # Strong attractive potential with repulsive forces, for A-B and C-D (the concatenation of the following two)
    starting_point_ppm = 1
    rmin = the_repulsive_rc * (3 / (1 + 2 * (the_repulsive_rc / the_sigma) ** 2)) ** (1 / 2)
    the_WF_potential_start_1_attract = np.array(
        [WF_potential_shift(shift=starting_point_ppm-1, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc,
                            r=i) + repulsive_barrier_epsilon
         for i in np.arange(starting_point_ppm, starting_point_ppm + (the_rc - 1), 1e-3)])
    the_WF_force_start_1_attract = np.array(
        [WF_force_shift(shift=starting_point_ppm-1, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange(starting_point_ppm, starting_point_ppm + (the_rc - 1), 1e-3)])

    the_WF_potential_start_1_repulse = np.array(
        [-WF_potential_shift(shift=the_rc - 1 - (rmin - 1), eps=repulsive_barrier_epsilon, sigma=the_sigma,
                             nu=1, mu=1, rc=the_rc, r=i)
         for i in
         np.arange(starting_point_ppm + (the_rc - 1), starting_point_ppm + (the_rc - 1) + the_rc - 1 - (rmin - 1),
                   1e-3)])
    the_WF_force_start_1_repulse = np.array(
        [-WF_force_shift(shift=the_rc - 1 - (rmin - 1), eps=repulsive_barrier_epsilon, sigma=the_sigma,
                         nu=1, mu=1, rc=the_rc, r=i)
         for i in
         np.arange(starting_point_ppm + (the_rc - 1), starting_point_ppm + (the_rc - 1) + the_rc - 1 - (rmin - 1),
                   1e-3)])
    PM_strong_potential = np.append(the_WF_potential_start_1_attract, the_WF_potential_start_1_repulse)
    PM_strong_force = np.append(the_WF_force_start_1_attract, the_WF_force_start_1_repulse)



    WF0 = hoomd.md.pair.Table(nlist=cell)
    # Reject big bodies
    WF0.params[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = dict(
        r_min=(diameter_particle - 1) + 0.01,
        U=WF_reject_particle_potential,
        F=WF_reject_particle_force)
    WF0.r_cut[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = (diameter_particle - 1) + 1.0
    ##############
    ##############
    PM_weak = hoomd.md.pair.Table(nlist=cell)
    PM_weak.params[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = dict(
        r_min=1,
        U=PM_weak_potential,
        F=PM_weak_force)
    PM_weak.r_cut[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = 0
    ##############
    ##############
    PM_strong = hoomd.md.pair.Table(nlist=cell)
    PM_strong.params[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = dict(
        r_min=1,
        U=PM_strong_potential,
        F=PM_strong_force)
    PM_strong.r_cut[(("C1", "C2", "A", "B"), ("C1", "C2", "A", "B"))] = 0
    ##############
    ##############
    ############################################
    ###########################################
    PM_weak.params[("C1", "A")] = dict(
        r_min=1,
        U=PM_weak_potential,
        F=PM_weak_force)
    PM_weak.r_cut[("C1", "A")] = diameter_particle + (the_rc - 1)
    PM_weak.params[("C2", "B")] = dict(
        r_min=1,
        U=PM_weak_potential,
        F=PM_weak_force)
    PM_weak.r_cut[("C2", "B")] = diameter_particle + (the_rc - 1)
    PM_strong.params[("C1", "C2")] = dict(
        r_min=1,
        U=PM_strong_potential,
        F=PM_strong_force)
    PM_strong.r_cut[("C1", "C2")] = starting_point_ppm + (the_rc - 1) + (the_repulsive_rc - 1) - (rmin - 1)
    PM_strong.params[("A", "B")] = dict(
        r_min=1,
        U=PM_strong_potential,
        F=PM_strong_force)
    PM_strong.r_cut[("A", "B")] = starting_point_ppm + (the_rc - 1) + (the_repulsive_rc - 1) - (rmin - 1)


    integrator.forces.append(WF0)
    integrator.forces.append(PM_weak)
    integrator.forces.append(PM_strong)
    nvt1 = hoomd.md.methods.Brownian(kT=the_kT, filter=hoomd.filter.All(), alpha=(1/diameter_particle)*drag_coeff)
    integrator.methods.append(nvt1)

    #############################################
    ####         Initialize simulation      ####
    #############################################
    if device == "cpu":
        cpu = hoomd.device.CPU()
        sim = hoomd.Simulation(device=cpu, seed=run)
    else:
        gpu = hoomd.device.GPU()
        sim = hoomd.Simulation(device=gpu, seed=run)
    sim.operations.integrator = integrator
    sim.create_state_from_gsd(filename=folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon))
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=the_kT)

    if saving_trajectory == True:
        gsd_writer = hoomd.write.GSD(filename=folder_to_store + '/MD' +
        '/traj_length_{}_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(length_box,np.round(run, 2),
        np.round(weak_epsilon,2),
        np.round(strong_epsilon,2),
        np.round(repulsive_barrier_epsilon, 2)),
        trigger=hoomd.trigger.Periodic(the_recoding_period),
        mode='wb')
        sim.operations.writers.append(gsd_writer)

    if stop_at_event == False:
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All())
        sim.operations.computes.append(thermodynamic_properties)
        sim.run(0)
        print("Initial potential: ", thermodynamic_properties.potential_energy)
        a_step = 0
        while a_step < nb_steps:
            sim.run(1)
            a_step +=1
        print(thermodynamic_properties.potential_energy)
        os.remove(folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon))

    else:
        a_step = 0
        substrate_free = 0
        broke = 0
        state3 = 0
        list_2d = [length_box, length_box]
        while a_step <= nb_steps and substrate_free == 0:
            a_step += 1
            sim.run(length_step)
            a_snapshot = sim.state.get_snapshot()
            snapshot_position = a_snapshot.particles.position
            posA = snapshot_position[2][:-1]
            posB = snapshot_position[3][:-1]
            posC1 = snapshot_position[0][:-1]
            posC2 = snapshot_position[1][:-1]
            if distance(posC1, posC2, list_2d) > diameter_particle + (the_rc - 1)+0.1:
                print("CATALYST BROKE at run:", run) # it will also be reported in "broke".
                broke = 1
                substrate_free = 1
                os.remove(
                    folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(
                        run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon))
                return broke, substrate_free, state3, the_dt * a_step * length_step

            if distance(posA, posB, list_2d) < diameter_particle + (the_rc - 1):
                if state3 ==0: # Check when they formed they bound (just as an additional information)
                    if distance(posA, posB, list_2d) < 1.05:
                        if distance(posA, posC1, list_2d) < 1.05:
                            if distance(posB, posC2, list_2d) < 1.05:
                                state3 = the_dt * a_step * length_step
                if distance(posA, posC1, list_2d) > number_beyond:
                    if distance(posB, posC2, list_2d) > number_beyond:
                        substrate_free = 1  # stop the simulation
                        os.remove(
                            folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon))
                        return broke, state3, the_dt * a_step * length_step
        os.remove(folder_to_store + '/MD' + '/init_run_{}_weak_epsilon_{}_strong_epsilon_{}_repulsiveEps_{}.gsd'.format(run, weak_epsilon, strong_epsilon, repulsive_barrier_epsilon))
        return broke, state3, the_dt * a_step * length_step


broke, state3, time =Model1(folder_to_store='/Users/yannsakref/Data/MOLECULAR_DYNAMICS/Energetic/Model1/test',
                                    length_box=5, diameter_particle=1,
                                    list_positions=[[0, 0, 0], [1.05, 0, 0], [0, 1.05, 0], [1.05, 1.05, 0]],
                                    list_orientation=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]],
                                    weak_epsilon=0.5, strong_epsilon=40, repulsive_barrier_epsilon=1, the_dt=1e-5, nb_steps=1e5, length_step = 100,
                                    the_recoding_period=100, run=0, drag_coeff=10, the_kT=1,
                                    stop_at_event=True, saving_trajectory=True, device = "cpu", the_rc = 1.1, the_repulsive_rc = 1.1, number_beyond = 1.1)
print(broke, state3, time)
