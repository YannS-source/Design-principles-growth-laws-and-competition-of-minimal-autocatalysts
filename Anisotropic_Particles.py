import sys
sys.path.append("/usr/local")
import gsd.hoomd
import hoomd
import numpy as np
import os
import cmath

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
# function to determine the position of patches.
def calc(count, diameter, center):
    x, y = center
    for i in range(count):
        r = cmath.rect(diameter / 2, (2 * cmath.pi) * (i / count))
        yield [round(x + r.real, 2), round(y + r.imag, 2), 0]


def Model2_v4(folder_to_store='/Users/yannsakref/Data/MOLECULAR_DYNAMICS', scheme_type="1A",
              length_box=5, diameter_patch=0.1,
              list_positions=[[0, 0, 0], [1.05, 0, 0], [0, 1.05, 0], [1.05, 1.05, 0]],
              list_orientation=[[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0],
                                [0.0, 0.0, 1.0, 0.0]], tilted_position=0,
              weak_epsilon=0.5, strong_epsilon=40, repulsive_barrier_epsilon=1, the_dt=1e-5,
              nb_steps=1e5, length_step=100,
              the_recoding_period=100, run=0, brownian=True, drag_coeff=10, the_kT=1,
              stop_at_event=False, saving_trajectory=False, device="cpu", the_rc=1.1,
              the_repulsive_rc=1.1):
    """
    Input:
    length_box: length of the square box, with periodic boundary conditions.
    diameter_patch: diameter of the patches on A, B, C, and D -- while diameter bodies are set to 1.
    list_position: list of the particle's position at initialization, [[A], [B], [C], [D]],
    where [i] means [x_i, y_i, z_i].
    list_orientation: list of the particle's position at start. Keep as such to ensure that the geometries are
     complementary.
    tilted_position: angle for the vertical patches. a default tilt of 0 means that the patches will be aligned, at
    90 degrees.
    weak_epsilon: interaction strength between the patches of A-C and B-D
    strong_epsilon: interaction strength between patches of A-B and C-D
    repulsive_barrier_epsilon: interaction barrier between A and B.
    the_dt: time step
    nb_steps: number of steps after which the simulation stops.
    length_step: length of the steps
    the_recording_period: in case saving_trajectory=True, the number of step between saving snapshots.
    run: random seed used.
    brownian: if true, Brownian motion, else, Langevin.
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
    num_bodies = 4 # A, B, C, and D (the small patches are not considered "body").
    L = length_box
    diameter_patch = diameter_patch
    diameter = 1
    snapshot = gsd.hoomd.Snapshot()
    snapshot.particles.N = num_bodies
    snapshot.particles.diameter = [diameter] * num_bodies
    snapshot.particles.position = list_positions
    snapshot.particles.moment_inertia = [[0, 0, 1]] * num_bodies
    snapshot.particles.types = ["C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A",
                                "W"] # C1 = C, C2 = D, A=A, B=B; the other particles are the patches
    snapshot.configuration.box = [L, L, 0, 0, 0, 0]
    snapshot.particles.orientation = list_orientation
    snapshot.particles.typeid = [0, 1, 2, 3]

    with gsd.hoomd.open(
            name=folder_to_store + "/MD" + '/init_length_{}_barrier_{}_weak_{}_tilt_{}_run_{}.gsd'.format(
                    length_box, repulsive_barrier_epsilon,
                    weak_epsilon, tilted_position,
                    run), mode='wb') as f:
        f.append(snapshot)

    ##############################################
    #### Precise the structure of the bodies  ####
    #############################################
    rigid = hoomd.md.constrain.Rigid()

    list_angles = [0, np.pi / 2 + tilted_position] # we are going to consider only two patches per rigid body
    list_of_positions_particles_C1 = [compute_position_based_on_angle(diameter, i) for i in list_angles]
    rigid.body['C1'] = {
        "constituent_types": ["C1-C2"] * 1 + ["C1-A"] * 1,
        "positions": list_of_positions_particles_C1,
        "orientations": [(1.0, 0.0, 0.0, 0.0)] * len(list_of_positions_particles_C1),
        "charges": [0.0] * len(list_of_positions_particles_C1),
        "diameters": [diameter_patch] * len(list_of_positions_particles_C1)}
    list_of_positions_particles_C2 = [compute_position_based_on_angle(diameter, i) for i in list_angles]
    rigid.body['C2'] = {
        "constituent_types": ["C2-C1"] * 1 + ["C2-B"] * 1,
        "positions": list_of_positions_particles_C2,
        "orientations": [(1.0, 0.0, 0.0, 0.0)] * len(list_of_positions_particles_C2),
        "charges": [0.0] * len(list_of_positions_particles_C2),
        "diameters": [diameter_patch] * len(list_of_positions_particles_C2)}
    list_angles = [0, (3 * np.pi / 2) - tilted_position]
    list_of_positions_particles_A = [compute_position_based_on_angle(diameter, i) for i in list_angles]
    rigid.body['A'] = {
        "constituent_types": ['A-B'] * 1 + ["A-C1"] * 1,
        "positions": list_of_positions_particles_A,
        "orientations": [(1.0, 0.0, 0.0, 0.0)] * len(list_of_positions_particles_A),
        "charges": [0.0] * len(list_of_positions_particles_A),
        "diameters": [diameter_patch] * len(list_of_positions_particles_A)}
    rigid.body['B'] = {
        "constituent_types": ['B-A'] * 1 + ["B-C2"] * 1,
        "positions": list_of_positions_particles_A,
        "orientations": [(1.0, 0.0, 0.0, 0.0)] * len(list_of_positions_particles_A),
        "charges": [0.0] * len(list_of_positions_particles_A),
        "diameters": [diameter_patch] * len(list_of_positions_particles_A)}

    #############################################
    ####         Initialized integrator      ####
    #############################################
    integrator = hoomd.md.Integrator(dt=the_dt, integrate_rotational_dof=True)
    integrator.rigid = rigid
    cell = hoomd.md.nlist.Cell(buffer=0, exclusions=['body'])

    # the same potentials as for the case of isotropic particles.
    def WF_potential_shift(shift=0, eps=1, sigma=1, nu=1, mu=1, rc=1.2, r=1):
        alpha = 2 * nu * (rc / sigma) ** (2 * mu) * ((1 + 2 * nu) / ((2 * nu) * ((rc / sigma) ** (2 * mu) - 1))) ** (
                2 * nu + 1)
        return eps * alpha * ((sigma / (r - shift)) ** (2 * mu) - 1) * ((rc / (r - shift)) ** (2 * mu) - 1) ** (2 * nu)
    def WF_force_shift(shift=0, eps=1, sigma=1, nu=1, mu=1, rc=1.2, r=1):
        alpha = 2 * nu * (rc / sigma) ** (2 * mu) * ((1 + 2 * nu) / ((2 * nu) * ((rc / sigma) ** (2 * mu) - 1))) ** (
                2 * nu + 1)
        CC = (rc / (r - shift)) ** (2 * mu)
        SS = (sigma / (r - shift)) ** (2 * mu)
        return (2 / (r - shift)) * eps * alpha * mu * (CC - 1) ** (2 * nu - 1) * (
                    2 * nu * (SS - 1) * CC + SS * (CC - 1))

    the_sigma = 1
    the_repulsive_epsilon = 1
    WF_reject_particle_potential = np.array(
        [WF_potential_shift(shift=0, eps=the_repulsive_epsilon, sigma=the_sigma, nu=1, mu=1, rc=1.1, r=i) for i in
         np.arange(0.01, 1.00, 0.01)])
    WF_reject_particle_force = np.array(
        [WF_force_shift(shift=0, eps=the_repulsive_epsilon, sigma=the_sigma, nu=1, mu=1, rc=1.1, r=i) for i in
         np.arange(0.01, 1.00, 0.01)])

    the_WF_potential_start_1_weak = np.array(
        [WF_potential_shift(shift=-0.9999, eps=weak_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange(1 - 0.9999, the_rc - 0.9999, 1e-3)])
    the_WF_force_start_1_weak = np.array(
        [WF_force_shift(shift=-0.9999, eps=weak_epsilon, sigma=the_sigma, nu=1, mu=1, rc=the_rc, r=i) for i in
         np.arange(1 - 0.9999, the_rc - 0.9999, 1e-3)])

    starting_point_ppm = 1
    rmin = the_repulsive_rc * (3 / (1 + 2 * (the_repulsive_rc / the_sigma) ** 2)) ** (1 / 2)
    the_WF_potential_start_1_attract = np.array(
        [WF_potential_shift(shift=starting_point_ppm - 1 - 0.9999, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1,
                            rc=the_rc,
                            r=i) + repulsive_barrier_epsilon
         for i in np.arange(starting_point_ppm - 0.9999, starting_point_ppm - 0.9999 + (the_rc - 1), 1e-3)])
    the_WF_force_start_1_attract = np.array(
        [WF_force_shift(shift=starting_point_ppm - 1 - 0.9999, eps=strong_epsilon, sigma=the_sigma, nu=1, mu=1,
                        rc=the_rc, r=i) for i in
         np.arange(starting_point_ppm - 0.9999, starting_point_ppm - 0.9999 + (the_rc - 1), 1e-3)])

    the_WF_potential_start_1_repulse = np.array(
        [-WF_potential_shift(shift=the_rc - 1 - (rmin - 1) - 0.9999, eps=repulsive_barrier_epsilon, sigma=the_sigma,
                             nu=1, mu=1, rc=the_rc, r=i)
         for i in
         np.arange(starting_point_ppm - 0.9999 + (the_rc - 1),
                   starting_point_ppm - 0.9999 + (the_rc - 1) + the_rc - 1 - (rmin - 1),
                   1e-3)])
    the_WF_force_start_1_repulse = np.array(
        [-WF_force_shift(shift=the_rc - 1 - (rmin - 1) - 0.9999, eps=repulsive_barrier_epsilon, sigma=the_sigma,
                         nu=1, mu=1, rc=the_rc, r=i)
         for i in
         np.arange(starting_point_ppm - 0.9999 + (the_rc - 1),
                   starting_point_ppm - 0.9999 + (the_rc - 1) + the_rc - 1 - (rmin - 1),
                   1e-3)])


    PP_strong_potential = np.append(the_WF_potential_start_1_attract, the_WF_potential_start_1_repulse)
    PP_strong_force = np.append(the_WF_force_start_1_attract, the_WF_force_start_1_repulse)


    WF0 = hoomd.md.pair.Table(nlist=cell)

    # Reject big bodies
    WF0.params[(("C1", "C2", "A", "B"),
                ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = dict(
        r_min=(1 - 1) + 0.01,
        U=WF_reject_particle_potential,
        F=WF_reject_particle_force)
    WF0.r_cut[(("C1", "C2", "A", "B"),
               ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = 1.0
    ####
    WF0.params[(("C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"),
                ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = dict(
        r_min=(1 - 1) + 0.01,
        U=WF_reject_particle_potential,
        F=WF_reject_particle_force)
    WF0.r_cut[(("C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"),
               ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = 0
    ##############
    PP_weak = hoomd.md.pair.Table(nlist=cell)
    PP_weak.params[(("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"), (
    "C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = dict(
        r_min=1,
        U=the_WF_potential_start_1_weak,
        F=the_WF_force_start_1_weak)
    PP_weak.r_cut[(("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"),
                   ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = 0
    ###
    PP_strong = hoomd.md.pair.Table(nlist=cell)
    PP_strong.params[(("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"), (
    "C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = dict(
        r_min=1,
        U=PP_strong_potential,
        F=PP_strong_force)
    PP_strong.r_cut[(("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"),
                     ("C1", "C2", "A", "B", "C1-C2", "C2-C1", "C1-A", "A-C1", "C2-B", "B-C2", "A-B", "B-A", "W"))] = 0
    ##############

    ## Weak
    PP_weak.params[("C1-A", "A-C1")] = dict(
        r_min=1 - 0.9999,
        U=the_WF_potential_start_1_weak,
        F=the_WF_force_start_1_weak)
    PP_weak.r_cut[("C1-A", "A-C1")] = the_rc - 0.9999
    PP_weak.params[("C2-B", "B-C2")] = dict(
        r_min=1 - 0.9999,
        U=the_WF_potential_start_1_weak,
        F=the_WF_force_start_1_weak)
    PP_weak.r_cut[("C2-B", "B-C2")] = the_rc - 0.9999
    ## Strong
    PP_strong.params[("C1-C2", "C2-C1")] = dict(r_min=1 - 0.9999, U=PP_strong_potential, F=PP_strong_force)
    PP_strong.r_cut[("C1-C2", "C2-C1")] = starting_point_ppm - 0.9999 + (the_rc - 1) + the_repulsive_rc - 1 - (rmin - 1)
    PP_strong.params[("A-B", "B-A")] = dict(r_min=1 - 0.9999, U=PP_strong_potential, F=PP_strong_force)
    PP_strong.r_cut[("A-B", "B-A")] = starting_point_ppm - 0.9999 + (the_rc - 1) + the_repulsive_rc - 1 - (rmin - 1)

    integrator.forces.append(WF0)
    integrator.forces.append(PP_weak)
    integrator.forces.append(PP_strong)

    if brownian:
        nvt = hoomd.md.methods.Brownian(kT=the_kT, filter=hoomd.filter.Rigid(("center", "free")), alpha=drag_coeff)
    else:
        nvt = hoomd.md.methods.Langevin(kT=the_kT, filter=hoomd.filter.Rigid(("center", "free")), alpha=drag_coeff)

    integrator.methods.append(nvt)

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
    sim.create_state_from_gsd(
        filename=folder_to_store + "/MD" + '/init_length_{}_barrier_{}_weak_{}_tilt_{}_run_{}.gsd'.format(
            length_box, repulsive_barrier_epsilon,
            weak_epsilon, tilted_position,
            run))
    rigid.create_bodies(sim.state)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=the_kT)
    if saving_trajectory == True:
        gsd_writer = hoomd.write.GSD(
            filename=folder_to_store + "/MD" + "/traj_length_{}_barrier_{}_weak_{}_tilt_{}_beta_{}_run_{}.gsd".format(
                length_box, repulsive_barrier_epsilon,
                weak_epsilon, tilted_position, 1, run),
            trigger=hoomd.trigger.Periodic(the_recoding_period),
            mode='wb')
        sim.operations.writers.append(gsd_writer)


    if stop_at_event == False:
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All())
        sim.operations.computes.append(thermodynamic_properties)
        sim.run(0)
        print("Initial potential: ", thermodynamic_properties.potential_energy)
        sim.run(nb_steps)
        print("Final potential: ", thermodynamic_properties.potential_energy)
    else:
            a_step = 0
            substrate_free = 0
            moment_in_state3 = 0
            broke = 0
            list_2d = np.array([length_box, length_box, 0])
            while a_step <= nb_steps and substrate_free == 0:
                a_step += 1
                sim.run(length_step)
                a_snapshot = sim.state.get_snapshot()
                snapshot_position = a_snapshot.particles.position
                snapshot_id = list(a_snapshot.particles.typeid)
                posC1C2 = snapshot_position[snapshot_id.index(4)]
                posC2C1 = snapshot_position[snapshot_id.index(5)]
                if distance(posC1C2, posC2C1, list_2d) > 1 + (the_rc - 1):
                    print("CATALYST BROKE at run:", run)
                    broke = 1
                    substrate_free = 1
                    os.remove(
                        folder_to_store + '/MD' + '/init_length_{}_barrier_{}_weak_{}_tilt_{}_run_{}.gsd'.format(
                            length_box, repulsive_barrier_epsilon,
                            weak_epsilon, tilted_position,
                            run))
                    print("Have removed")
                    return broke, moment_in_state3, the_dt * a_step * length_step

                posAB = snapshot_position[snapshot_id.index(10)]
                posBA = snapshot_position[snapshot_id.index(11)]
                if distance(posAB, posBA, list_2d) < the_rc - diameter:  # AB and BA interacts
                    posAC1 = snapshot_position[snapshot_id.index(6)]
                    posC1A = snapshot_position[snapshot_id.index(7)]
                    if distance(posAC1, posC1A, list_2d) > the_rc - diameter:
                        posBC2 = snapshot_position[snapshot_id.index(8)]
                        posC2B = snapshot_position[snapshot_id.index(9)]
                        if distance(posBC2, posC2B, list_2d) > the_rc - diameter:
                            os.remove(
                                folder_to_store + "/MD" + '/init_length_{}_barrier_{}_weak_{}_tilt_{}_run_{}.gsd'.format(
                                    length_box, repulsive_barrier_epsilon,
                                    weak_epsilon, tilted_position,
                                    run))
                            print("A time:", the_dt * a_step * length_step)
                            print("Corresponding step", a_step)
                            print("That is, ", the_dt * a_step * length_step / (the_dt * length_step))
                            return broke, moment_in_state3, the_dt * a_step * length_step


    os.remove(
        folder_to_store + "/MD" + '/init_length_{}_barrier_{}_weak_{}_tilt_{}_run_{}.gsd'.format(length_box,
                                                                                                         repulsive_barrier_epsilon,
                                                                                                         weak_epsilon,
                                                                                                         tilted_position,
                                                                                                         run))
    return the_dt * a_step * length_step