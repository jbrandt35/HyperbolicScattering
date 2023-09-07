from HyperbolicTools import *
import time

start_time = time.time()

# setting parameters of simulation, set in Setup.py
#'''
perturber_a = PERTSEPARATION
inner_a = INNERSEMIMAJ
#'''

migration_type = "MIGRATION_TYPE"
migration_constant = MIGRATION_CONST
gr_type = "GR_TYPE"
randomize_M = RANDOMIZE_M
randomize_inc = RANDOMIZE_INC
inc_std = INCLINATION_STD
randomize_ecc = RANDOMIZE_ECC
randomize_mass = RANDOMIZE_MASS

m0 = SMBH_MASS
m1 = INNER_MASS
m2 = PERTURBER_MASS

m1_range = [INNER_MASS_LOWER_BOUND, INNER_MASS_UPPER_BOUND]
m2_range = [PERTURBER_MASS_LOWER_BOUND, PERTURBER_MASS_UPPER_BOUND]

if randomize_mass:
    m1 = (m1_range[1] - m1_range[0])*np.random.rand() + m1_range[0]
    m2 = (m2_range[1] - m2_range[0])*np.random.rand() + m2_range[0]

#####################################################################################

sim = create_simulation() # creating simulation, adding particles
inner_M, perturber_M, inner_inc, perturber_inc, inner_e, perturber_e = populate_simulation(sim, perturber_a = perturber_a, randomize_M = randomize_M, randomize_inc = randomize_inc, inc_std=inc_std, randomize_ecc = randomize_ecc, SMBH_m=m0, inner_m=m1, perturber_m=m2)


inner_period, perturber_period = get_inner_period(sim), get_perturber_period(sim) # getting periods, determine simulation time, timestep
sim.dt = 0.05 * inner_period
sim_time = 1e5 * inner_period
SMBH_Rg = sim.G * m0 / c ** 2
inner_a *= SMBH_Rg # initial semimaj axis, needed for trap radius, 9.870628716888302 AU for 1000 Rg of SMBH initially
distance_limit = 100*inner_a # simulation aborts if particles become farther than this from SMBH
inner_Rg = sim.G * m1 / c ** 2
perturber_Rg = sim.G * m2 / c ** 2

def dragForce(reb_sim, rebx_force, particles, N): # drag force from Li & Lai
    sim = reb_sim.contents
    dragforce = rebx_force.contents
    vxyz = particles[2].vxyz
    r_0 = particles[0].xyz
    r_pert = particles[2].xyz
    tau_drag = dragforce.params["c"]
    ax, ay, az = calc_drag_force(tau_drag*inner_period, vxyz, sim.G, particles[0].m, r_0, r_pert)
    particles[2].ax += ax
    particles[2].ay += ay
    particles[2].az += az

def trapForce(reb_sim, rebx_force, particles, N): # trap force from Li & Lai, have not used in testing for a while
    sim = reb_sim.contents
    trapforce = rebx_force.contents
    r_SMBH = particles[0].xyz
    r_pert = particles[2].xyz
    tau_trap = trapforce.params["c"]
    ax, ay, az = calc_trap_force(tau_trap*inner_period, sim.G, particles[0].m, inner_a, r_pert, r_SMBH)
    particles[2].ax += ax
    particles[2].ay += ay
    particles[2].az += az
    
rebx = reboundx.Extras(sim)

# old way to include rebx effects, new way is done through setup file
'''
if include_drag_force:
    myforce = rebx.create_force("drag")
    myforce.force_type = "vel"
    myforce.update_accelerations = dragForce
    rebx.add_force(myforce)
    myforce.params["c"] = tau_drag

if include_trap_force:
    myforce = rebx.create_force("trap")
    myforce.force_type = "vel"
    myforce.update_accelerations = trapForce
    rebx.add_force(myforce)
    myforce.params["c"] = tau_trap

if include_rebound_migration:
    mig = rebx.load_force("modify_orbits_forces")
    rebx.add_force(mig)
    sim.particles['perturber'].params['tau_a'] = -tau_mig
    print('Migration force initialized')

if include_gr_old:
    
    gr_full = rebx.load_force("gr_full")
    rebx.add_force(gr_full)
    gr_full.params["c"] = c
    
    gr_radiation = rebx.load_force("gr_radiation") # gravitational wave radiation, not included right now as migration being tested
    rebx.add_force(gr_radiation)
    gr_radiation.params["c"] = c
    gr_radiation.params["gr_rad_part1"] = 1
    gr_radiation.params["gr_rad_part2"] = 2
#'''

if migration_type == "drag":
    tau_drag = migration_constant
    myforce = rebx.create_force("drag")
    myforce.force_type = "vel"
    myforce.update_accelerations = dragForce
    rebx.add_force(myforce)
    myforce.params["c"] = tau_drag
    
if migration_type == "trap":
    tau_trap = migration_constant
    sim = reb_sim.contents
    trapforce = rebx_force.contents
    r_SMBH = particles[0].xyz
    r_pert = particles[2].xyz
    tau_trap = trapforce.params["c"]
    ax, ay, az = calc_trap_force(tau_trap*inner_period, sim.G, particles[0].m, inner_a, r_pert, r_SMBH)
    particles[2].ax += ax
    particles[2].ay += ay
    particles[2].az += az

# different post-Newtonian effects

if '1PN' in gr_type: # 1PN from REBOUND
    gr_full = rebx.load_force("gr_full")
    rebx.add_force(gr_full)
    gr_full.params["c"] = c

if '2.5PN' in gr_type: # from Hareesh, just 2.5PN
    gr_radiation = rebx.load_force("gr_radiation")
    rebx.add_force(gr_radiation)
    gr_radiation.params["c"] = c
    gr_radiation.params["gr_rad_part1"] = 1
    gr_radiation.params["gr_rad_part2"] = 2

'''
if gr_type == "new": # 1PN+2.5PN, only active when particles within 100Rg
    gr_radiation_full = rebx.load_force("gr_radiation_full")
    rebx.add_force(gr_radiation_full)
    gr_radiation_full.params["c"] = c
    gr_radiation_full.params["gr_rad_part1"] = 1
    gr_radiation_full.params["gr_rad_part2"] = 2
#'''


initialize_data_collection()
#####################################################################################

outcome_record["inner_M"] = inner_M
outcome_record["perturber_M"] = perturber_M
outcome_record["inner_inc"] = inner_inc
outcome_record["perturber_inc"] = perturber_inc
outcome_record["inner_e"] = inner_e
outcome_record["perturber_e"] = perturber_e
outcome_record["inner_a"] = inner_a
outcome_record["perturber_a"] = perturber_a
outcome_record['SMBH_m'] = m0
outcome_record['inner_m'] = m1
outcome_record['perturber_m'] = m2

dump_record(outcome_record)

min_dist_dict = {'Time (%)': [], 'Time (yr)': [], 'Step': [], 'Separation (AU)': [], 'Separation (Rg)': [], 'True Anomaly': [],'Inner x': [], 'Inner y': [], 'Inner z': [], 'Inner vx': [], 'Inner_vy': [], 'Inner_vz': [], 'Perturber x': [], 'Perturber y': [], 'Perturber z': [], 'Perturber vx': [], 'Perturber vy': [], 'Perturber vz': [],}
rolling_df = pd.DataFrame(min_dist_dict)
rolling_df_size = 5000

global total_time_steps_completed
'''
while total_time_steps_completed < rolling_df_size:

    time_percent = sim.t/sim_time
    d = dist_between(sim.particles["perturber"], sim.particles["inner"])

    rolling_df.loc[len(rolling_df.index)] = [time_percent, sim.t, total_time_steps_completed, d, d/inner_Rg, sim.particles["inner"].xyz[0], sim.particles["inner"].xyz[1], sim.particles["inner"].xyz[2], sim.particles["inner"].vxyz[0], sim.particles["inner"].vxyz[1], sim.particles["inner"].vxyz[2], sim.particles["perturber"].xyz[0], sim.particles["perturber"].xyz[1], sim.particles["perturber"].xyz[2], sim.particles["perturber"].vxyz[0], sim.particles["perturber"].vxyz[1], sim.particles["perturber"].vxyz[2],]

    total_time_steps_completed += 1
    sim.step()
#'''

#outcome_record["Minimum Separation"] = d
close_encounter_threshold = 50
outcome_record["Minimum Separation (AU)"] = close_encounter_threshold*(inner_Rg + perturber_Rg)
outcome_record["Minimum Separation (Rg)"] = close_encounter_threshold
new_event = False
event_index = 0
event_df = pd.DataFrame()

#'''
while sim.t <= sim_time: # should be 1e5 periods

    time_percent = sim.t/sim_time
    d = dist_between(sim.particles["perturber"], sim.particles["inner"])
       
    if not outcome_record["Formed Binary"]:
        if is_bound(sim.particles["perturber"], sim.particles["inner"]):
            outcome_record["Formed Binary"] = True
    
    if total_time_steps_completed % 1000 == 0:
        save_data(sim)

    if new_event: # saving encounter in progress
        append_data(event_df, time_percent, sim.t, total_time_steps_completed, d, d/(inner_Rg + perturber_Rg), sim)
        if d < min_d:
            min_d = d
        elif d > close_encounter_threshold*(inner_Rg + perturber_Rg): # ending close encounter
            new_event = False
            event_df.to_csv('close_encounters/' + str(event_index) + '.csv')
            if min_d < outcome_record['Minimum Separation (AU)']:
                outcome_record['Minimum Separation (AU)'] = min_d

    elif d < close_encounter_threshold*(inner_Rg + perturber_Rg): # trigger new close encounter
        new_event = True
        min_d = d
        event_index += 1
        event_df = rolling_df.copy()

    #check_and_assign_minimums(sim, sim_time, outcome_record, df, total_time_steps_completed)
    final_time = time.time()
    elapsed = final_time-start_time
    outcome_record['Time'] = elapsed
    if d < outcome_record['Minimum Separation (AU)']:
        outcome_record['Minimum Separation (AU)'] = d
    check_bound(sim, outcome_record, distance_limit, start_time) # checking exit conditions
    check_for_collisions(sim, outcome_record, new_event, event_df, event_index, start_time)

    total_time_steps_completed += 1
    sim.step()
#'''

# old way to save encounters
'''
while sim.t <= sim_time: # should be 1e5
#while total_time_steps_completed < 55:

    time_percent = sim.t/sim_time
    d = dist_between(sim.particles["perturber"], sim.particles["inner"])

    append_data(rolling_df, time_percent, sim.t, total_time_steps_completed, d, d/inner_Rg, sim)
    rolling_df = rolling_df.drop([rolling_df.index[0]])
       
    if not outcome_record["Formed Binary"]:
        if is_bound(sim.particles["perturber"], sim.particles["inner"]):
            outcome_record["Formed Binary"] = True
    
    if total_time_steps_completed % 1000 == 0:
        save_data(sim)

    if min_reached:
        if event_step < rolling_df_size:
            append_data(event_df, time_percent, sim.t, total_time_steps_completed, d, d/inner_Rg, sim)
            event_step += 1
        else:
            min_reached, new_event = False, False
            event_df.to_csv('close_encounters/' + str(event_index) + '.csv')

    elif new_event:
        if d < last_d:
            last_d = d
        elif d > last_d:
            outcome_record['Minimum Separation (AU)'] = d
            outcome_record['Minimum Separation (Rg)'] = d/inner_Rg
            outcome_record['Minimum Separation Step'] = total_time_steps_completed
            event_df = rolling_df.copy()
            min_reached = True
            event_step = 0

    elif d < outcome_record['Minimum Separation (AU)']:
        new_event = True
        last_d = d
        event_index += 1

    #check_and_assign_minimums(sim, sim_time, outcome_record, df, total_time_steps_completed)

    check_bound(sim, outcome_record, distance_limit)
    check_for_collisions(sim, outcome_record, new_event, event_df, event_index)

    total_time_steps_completed += 1
    sim.step()

#'''
#####################################################################################

outcome_record["Result"] = "Ended with all particles within limit to SMBH" # if ends stable
if new_event:
    event_df.to_csv('close_encounters/' + str(event_index) + '.csv')
    if min_d < outcome_record['Minimum Separation (AU)']:
        outcome_record['Minimum Separation (AU)'] = min_d

save_final_data(outcome_record, sim, start_time)
dump_record(outcome_record)


