import json
import rebound
import numpy as np
from numpy.linalg import norm as mag
from astropy import constants
import pandas as pd
import os
import reboundx
import h5py
import matplotlib.pyplot as plt
import time

data_objects = ["inner", "SMBH", "perturber"]

#m0, m1, m2 = 10 ** 6, 20, 20 # solar mass
'''
m0 = SMBH_MASS
m1 = INNER_MASS
m2 = PERTURBER_MASS
#'''
#m0, m1_a, m1_b, m2 = 10 ** 6, 0.000000001, 0.000000001, 20
c = constants.c.to("AU/yr").value
au_to_lyr = 0.0000158125

#Returns a rebound simulation object with the appropriate units
def create_simulation():
    sim = rebound.Simulation()
    #sim.integrator = "whfast"
    sim.units = ["Msun", "yr", "AU"]
    return sim

sim = create_simulation()
G = sim.G
#inner_Rg = G * m1 / c ** 2
#print(inner_Rg)

# Returns the velocity of the COM of the binary given its orbital elements. Note: SMBH needs to be in sim.particles
def get_binary_COM_data(m_SMBH, m_binary, a, e = 0, M = 0, inc = 0):

    #Create a simulation just for getting these values
    virtual_sim = create_simulation()

    #Put in the SMBH
    virtual_sim.add(m = m_SMBH, hash = "SMBH")

    #Add a virtual particle to represent the COM of the binary
    virtual_sim.add(primary = virtual_sim.particles["SMBH"], m = m_binary, a = a, e = e, M = M, inc = inc, hash = "binary_COM")

    #Get its position
    position = np.array(virtual_sim.particles["binary_COM"].xyz)

    #Get its velocity
    velocity = np.array(virtual_sim.particles["binary_COM"].vxyz)

    #Get its Hill Radius
    R_hill = virtual_sim.particles["binary_COM"].rhill

    return position, velocity, R_hill


def get_BBH_data(m1, m2, a, e = 0, M = 0, inc = 0):

    # Create a simulation just for getting these values
    virtual_sim = create_simulation()

    #Add m1 as the primary object
    virtual_sim.add(m = m1, hash = "BH_1")

    #Add m2 around it
    virtual_sim.add(primary = virtual_sim.particles["BH_1"], m = m2, a = a, e = e, M = M, inc = inc, hash = "BH_2")

    #Switch to COM frame
    virtual_sim.move_to_com()

    #Get positions
    BH_1_position = np.array(virtual_sim.particles["BH_1"].xyz)
    BH_2_position = np.array(virtual_sim.particles["BH_2"].xyz)

    #Get velocities
    BH_1_velocity = np.array(virtual_sim.particles["BH_1"].vxyz)
    BH_2_velocity = np.array(virtual_sim.particles["BH_2"].vxyz)

    BH_1_data = BH_1_position, BH_1_velocity
    BH_2_data = BH_2_position, BH_2_velocity

    return BH_1_data, BH_2_data


#Takes in sim object and orbital elements, populates the simulation with the objects accordingly
#SMBH_a in units of Gravitational Radii, binary_separation in units of m1 Hill radii, perturber_a in units of m1+m2 Hill radii
#Returns w - for calculating mergers between m1 and m2
#Each particle is hashed with a name useful for other functions - i.e. see period calculation functions
'''
def populate_simulation(sim, binary_separation = 0.1, binary_eccentricity = 0, binary_M = 0, binary_inc = 0, SMBH_a = 1000, SMBH_eccentricity = 0, SMBH_M = 0, SMBH_inc = 0, perturber_a = 20, perturber_e = 0, perturber_M = 0, perturber_inc = 0, randomize_M = False, randomize_inc = False):

    #Calculate Gravitational Radius and convert units of SMBH_a
    Rg = sim.G * m0 / c ** 2
    SMBH_a *= Rg

    #Calculate m1+m2 Hill radius and find the semi-major axis of the perturber
    mass_factor = np.power((m1_a+m1_b+m2)/(3*m0), 1/3)
    perturber_factor = (1 + (perturber_a/2) * mass_factor) / (1 - (perturber_a/2) * mass_factor)
    perturber_a = SMBH_a * perturber_factor

    #If randomized mean anamolies are wanted, generate three numbers 0-2pi and set mean anamolies
    if randomize_M:
        SMBH_M, binary_M, perturber_M = 2 * np.pi * np.random.rand(3)
    
    mu, sigma = 0, 0.05
    if randomize_inc:
        perturber_inc = abs((np.random.normal(mu, sigma, 1))[0])
        perturber_inc *= np.pi/180

    #Add the SMBH to the sim
    sim.add(m = m0, hash = "SMBH")

    #Get position, velocity, and Hill radius of binary's COM
    binary_COM_position, binary_COM_velocity, m1_R_hill = get_binary_COM_data(m0, m1_a + m1_b, SMBH_a, e = SMBH_eccentricity, M = SMBH_M, inc = SMBH_inc)

    #Change units of binary_separation
    binary_separation *= m1_R_hill

    #Get data of BHs in binary (in binary COM F.O.R.)
    BH_1_data, BH_2_data = get_BBH_data(m1_a, m1_b, binary_separation, e = binary_eccentricity, M = binary_M, inc = binary_inc)

    #Unpack BBH data
    BH_1_position, BH_1_velocity = BH_1_data
    BH_2_position, BH_2_velocity = BH_2_data

    #Create particles for BBHs
    BH_1 = rebound.Particle(m = m1_a, hash = "BBH_1")
    BH_1.xyz = BH_1_position + binary_COM_position
    BH_1.vxyz = BH_1_velocity + binary_COM_velocity

    BH_2 = rebound.Particle(m = m1_b, hash = "BBH_2")
    BH_2.xyz = BH_2_position + binary_COM_position
    BH_2.vxyz = BH_2_velocity + binary_COM_velocity

    #Add BBHs into sim
    sim.add(BH_1)
    sim.add(BH_2)

    #Add perturber into sim
    sim.add(primary = sim.particles["SMBH"], m = m2, a = perturber_a, e = perturber_e, M = perturber_M, inc = perturber_inc, hash = "perturber")

    #Move to COM frame of simulation
    sim.move_to_com()

    return [np.sqrt(sim.G * m0 / SMBH_a) - np.sqrt(sim.G * m0 / perturber_a) + np.sqrt(
        sim.G * (m1_a + m1_b) / (binary_separation ** 3)) * binary_separation * (m1_b/(m1_a+m1_b)), SMBH_M, binary_M, perturber_M, SMBH_inc, binary_inc, perturber_inc]
'''

def populate_simulation(sim, inner_a = 1000, inner_e = 0, inner_M = 0, inner_inc = 0, perturber_a = 5.0, perturber_e = 0, perturber_M = 0, perturber_inc = 0, randomize_M = False, randomize_inc = False, inc_std = 0.05, randomize_ecc = False, SMBH_m=1e6, inner_m=20, perturber_m=20):

    #Calculate Gravitational Radius and convert units of SMBH_a
    Rg = sim.G * SMBH_m / c ** 2
    inner_a *= Rg

    #Calculate m1+m2 Hill radius and find the semi-major axis of the perturber
    mass_factor = np.power((inner_m+perturber_m)/(3*SMBH_m), 1/3)
    perturber_factor = (1 + (perturber_a/2) * mass_factor) / (1 - (perturber_a/2) * mass_factor)
    perturber_a = inner_a * perturber_factor

    #If randomized mean anamolies are wanted, generate three numbers 0-2pi and set mean anamolies
    if randomize_M:
        inner_M, perturber_M = 2 * np.pi * np.random.rand(2)
    
    
    if randomize_inc:
        mu, sigma = 0, inc_std

        perturber_inc = abs((np.random.normal(mu, sigma, 1))[0])
        perturber_inc *= np.pi/180

        inner_inc = abs((np.random.normal(mu, sigma, 1))[0])
        inner_inc *= np.pi/180

    if randomize_ecc:
        mu, sigma = 0.05, 0.02
        
        perturber_e_positive = False
        while perturber_e_positive == False:
            perturber_e = (np.random.normal(mu, sigma, 1))[0]
            if perturber_e > 0:
                perturber_e_positive = True

        inner_e_positive = False
        while inner_e_positive == False:
            inner_e = (np.random.normal(mu, sigma, 1))[0]
            if inner_e > 0:
                inner_e_positive = True


    #Add the SMBH to the sim
    sim.add(m = SMBH_m, hash = "SMBH")

    #Add inner stellar mass blackhole
    sim.add(primary = sim.particles["SMBH"], m = inner_m, a = inner_a, e = inner_e, M = inner_M, inc = inner_inc, hash = "inner")

    #Add perturber into sim
    sim.add(primary = sim.particles["SMBH"], m = perturber_m, a = perturber_a, e = perturber_e, M = perturber_M, inc = perturber_inc, hash = "perturber")

    #Move to COM frame of simulation
    sim.move_to_com()

    return [inner_M, perturber_M, inner_inc, perturber_inc, inner_e, perturber_e]


def calc_drag_force(tau, vxyz, G, M, r_0, r_pert):
    vx, vy, vz = vxyz
    r = [r_pert[0] - r_0[0], r_pert[1] - r_0[1], r_pert[2] - r_0[2]]
    x, y, z = r
    d = np.sqrt(x**2 + y**2 + z**2)
    beta = np.sqrt(G*M/(d**5))
    ax = (-1/tau)*(vx + y*beta)
    ay = (-1/tau)*(vy - x*beta)
    az = (-1/tau)*(vz)
    return ax, ay, az

def calc_trap_force(tau, G, M, a_bin, r_pert, r_SMBH):
    omega = np.sqrt(G*M/(a_bin**3))
    r = [r_pert[0] - r_SMBH[0], r_pert[1] - r_SMBH[1], r_pert[2] - r_SMBH[2]]
    x, y, z = r
    d = np.sqrt(x**2 + y**2 + z**2)
    F_mag = -omega*(d - a_bin)/tau
    ax = F_mag*(-y/d)
    ay = F_mag*(x/d)
    az = 0
    return ax, ay, az

def get_binary_period(sim):
    orbit = sim.particles["BBH_1"].calculate_orbit(primary = sim.particles["BBH_2"])
    return orbit.P


def get_inner_period(sim):
    orbit = sim.particles["inner"].calculate_orbit(primary = sim.particles["SMBH"])
    return orbit.P


def get_perturber_period(sim):
    orbit = sim.particles["perturber"].calculate_orbit(primary = sim.particles["SMBH"])
    return orbit.P


def is_bound(primary, secondary): # ecc check
    orbit = secondary.calculate_orbit(primary = primary)
    return bool(orbit.e < 1)

def system_ejected(primary, secondary, limit): # distance check
    d = dist_between(primary, secondary)
    return bool(d > limit)

def t_GW(sim, particle1, particle2): # gravitational wave timescale, Peters 1963

    orbit = particle1.calculate_orbit(primary = particle2)

    m1 = particle1.m
    m2 = particle2.m
    a = orbit.a
    e = orbit.e

    first_term = -64/5
    second_term = (sim.G ** 3) * (m1 + m2) * m1 * m2 / ((c**5) * (a**3) * ((1-e**2)**(7/2)))
    third_term = 1 + ((73/24)*(e**2)) + ((37/94)*(e**4))

    da_dt = first_term * second_term * third_term

    t_gw = np.abs(a/da_dt)
    return t_gw


def check_for_collisions(sim, record, new_event, df, event_index, start_time): # event horizon check

    ##############Check Binary Orbits###############

    BH_1 = sim.particles["inner"]
    perturber = sim.particles["perturber"]

    sum_of_event_horizons = 2 * sim.G * (BH_1.m + perturber.m) / (c ** 2)

    orbit = BH_1.calculate_orbit(primary = perturber)
    distance, period = orbit.d, orbit.P

    if distance < sum_of_event_horizons:
        record["Result"] = f"Collision Encountered: The distance between black holes was {distance} AU when the sum of event horizon radii was {sum_of_event_horizons} AU."
        save_final_data(record, sim, start_time)
        dump_record(record)
        if new_event:
            df.to_csv('close_encounters/' + str(event_index) + '.csv')
            outcome_record['Minimum Separation (AU)'] = distance
        raise CollisionException(record["Result"])


def check_bound(sim, record, limit, start_time):
    '''
    if not is_bound(sim.particles["BBH_1"], sim.particles["BBH_2"]):
        record["Result"] = f"Unbound: Binary eccentricity reached {sim.particles['BBH_2'].calculate_orbit(primary = sim.particles['BBH_1']).e}"
        dump_record(record)
        raise UnboundException(record["Result"])
    '''

    SMBH, BH_1, perturber = sim.particles
    if system_ejected(SMBH, BH_1, limit):
        record["Result"] = f"Ejected: Inner black hole left the system, ecc: {BH_1.calculate_orbit(primary = SMBH).e}"
        save_final_data(record, sim, start_time)
        dump_record(record)
        raise UnboundException(record["Result"])

    if system_ejected(SMBH, perturber, limit):
        record["Result"] = f"Ejected: Outer black hole left the system, ecc: {perturber.calculate_orbit(primary = SMBH).e}"
        save_final_data(record, sim, start_time)
        dump_record(record)
        raise UnboundException(record["Result"])


def dist_between(particle_1, particle_2):
    particle_1_position, particle_2_position = np.array(particle_1.xyz), np.array(particle_2.xyz)
    return mag(particle_1_position - particle_2_position)

def get_GW_power(m1, m2, r_p):
    mu = (m1*m2)/(m1+m2)
    m12 = m1 + m2
    numerator = 85*np.pi*(G**(7/2))*(mu**2)*(m12**(5/2))
    denominator = 12*(np.sqrt(2))*(c**5)*(r_p**(7/2))
    E_GW = numerator/denominator
    return E_GW

'''
def get_perturber_binary_separation(binary_COM, perturber):

    virtual_sim = create_simulation()
    virtual_sim.add(m = m0)
    virtual_sim.add(binary_COM)

    a_binary_COM = virtual_sim.particles[1].calculate_orbit(primary = virtual_sim.particles[0]).a
    a_perturber = perturber.calculate_orbit().a

    m = a_perturber - a_binary_COM

    del virtual_sim

    mass_factor = np.power(3*m0 / (m1_a + m1_b + m2), 1/3)
    a_factor = (a_perturber - a_binary_COM) / (a_perturber + a_binary_COM)

    n = 2 * mass_factor * a_factor

    return n
'''

'''
def check_and_assign_minimums(sim, sim_time, record, df, timestep):
    # Get binary COM as a rebound.Particle
    
    binary_barycenter = sim.calculate_com(first = 1, last = 3)

    distance_perturber_to_binary_COM = dist_between(sim.particles["perturber"], binary_barycenter)
    distance_between_BBHs = sim.particles["BBH_2"].calculate_orbit(primary = sim.particles["BBH_1"]).d

    t_gw = t_GW(sim, sim.particles["BBH_1"], sim.particles["BBH_2"])
    relative_tGW = t_gw / get_binary_period(sim)

    if distance_between_BBHs < record["Minimum Distance Between BBHs"]:
        record["Minimum Distance Between BBHs"] = distance_between_BBHs
    if distance_perturber_to_binary_COM < record["Minimum Distance Between Binary COM and Perturber"]:
        record["Minimum Distance Between Binary COM and Perturber"] = distance_perturber_to_binary_COM
    if t_gw < record["Minimum t_GW"]:
        record["Minimum t_GW"] = t_gw
    if relative_tGW < record["Minimum relative t_GW"]:
        record["Minimum relative t_GW"] = relative_tGW
'''

def check_and_assign_minimums(sim, sim_time, record, df, timestep): # dont really use this function anymore

    time = sim.t/(1e5 * get_inner_period(sim))
    check_distance = dist_between(sim.particles["perturber"], sim.particles["inner"])
    if check_distance < record["Minimum Separation"]:
        record["Minimum Separation"] = check_distance
        record["Minimum Separation Time"] = time
        record["Inner xyz at Min Sep"] = sim.particles["inner"].xyz
        record["Inner vxyz at Min Sep"] = sim.particles["inner"].vxyz
        record["Perturber xyz at Min Sep"] = sim.particles["perturber"].xyz
        record["Perturber vxyz at Min Sep"] = sim.particles["perturber"].vxyz

        df.loc[len(df.index)] = [time, timestep, check_distance, sim.particles["inner"].xyz[0], sim.particles["inner"].xyz[1], sim.particles["inner"].xyz[2], sim.particles["inner"].vxyz[0], sim.particles["inner"].vxyz[1], sim.particles["inner"].vxyz[2], sim.particles["perturber"].xyz[0], sim.particles["perturber"].xyz[1], sim.particles["perturber"].xyz[2], sim.particles["perturber"].vxyz[0], sim.particles["perturber"].vxyz[1], sim.particles["perturber"].vxyz[2],]

    if is_bound(sim.particles["perturber"], sim.particles["inner"]):
        t_gw = t_GW(sim, sim.particles["perturber"], sim.particles["inner"])
        if t_gw < record["Minimum t_GW"]:
            record["Minimum t_GW"] = t_gw
            record["Minimum t_GW Time"] = time
            record["Inner xyz at Min t_GW"] = sim.particles["inner"].xyz
            record["Inner vxyz at Min t_GW"] = sim.particles["inner"].vxyz
            record["Perturber xyz at Min t_GW"] = sim.particles["perturber"].xyz
            record["Perturber vxyz at Min t_GW"] = sim.particles["perturber"].vxyz


def append_data(df, time_percent, time, step, d, d_Rg, sim): # ways to save data conveniently during close encounter events
    #df.loc[len(df.index)+df.index[0]] = [time_percent, time, step, d, d_Rg, sim.particles["inner"].xyz[0], sim.particles["inner"].xyz[1], sim.particles["inner"].xyz[2], sim.particles["inner"].vxyz[0], sim.particles["inner"].vxyz[1], sim.particles["inner"].vxyz[2], sim.particles["perturber"].xyz[0], sim.particles["perturber"].xyz[1], sim.particles["perturber"].xyz[2], sim.particles["perturber"].vxyz[0], sim.particles["perturber"].vxyz[1], sim.particles["perturber"].vxyz[2],]
    df.loc[len(df.index)] = [time_percent, time, step, d, d_Rg, (sim.particles["inner"].calculate_orbit(primary=sim.particles['perturber'])).f, sim.particles["inner"].xyz[0], sim.particles["inner"].xyz[1], sim.particles["inner"].xyz[2], sim.particles["inner"].vxyz[0], sim.particles["inner"].vxyz[1], sim.particles["inner"].vxyz[2], sim.particles["perturber"].xyz[0], sim.particles["perturber"].xyz[1], sim.particles["perturber"].xyz[2], sim.particles["perturber"].vxyz[0], sim.particles["perturber"].vxyz[1], sim.particles["perturber"].vxyz[2],]

def append_data_alt(df, time, step, d, sim):
    df.loc[len(df.index)] = [time, step, d, sim.particles["inner"].xyz[0], sim.particles["inner"].xyz[1], sim.particles["inner"].xyz[2], sim.particles["inner"].vxyz[0], sim.particles["inner"].vxyz[1], sim.particles["inner"].vxyz[2], sim.particles["perturber"].xyz[0], sim.particles["perturber"].xyz[1], sim.particles["perturber"].xyz[2], sim.particles["perturber"].vxyz[0], sim.particles["perturber"].vxyz[1], sim.particles["perturber"].vxyz[2],]

def initialize_data_collection():
    os.system("rm -r -f result")
    os.system("mkdir result")
    os.system("rm -r -f close_encounters")
    os.system("mkdir close_encounters")
    print('initialized data collection')


def save_to_frame(df, data): # data saving stuff
    df.loc[len(df.index)] = data


def save_data_to_buffer(sim):

    save_to_frame(buffers["Misc"], [sim.t, sim.dt, dist_between(sim.particles["inner"], sim.particles["perturber"]), sim.particles["perturber"].calculate_orbit(primary = sim.particles["SMBH"]).a, sim.particles["inner"].calculate_orbit(primary = sim.particles["SMBH"]).a,])

    for particle in data_objects:
        save_to_frame(buffers[particle], sim.particles[particle].xyz)


def save_data_to_disk():
    with pd.HDFStore("result/data.h5") as data_file:
        data_file.append("Misc", buffers["Misc"])
        for particle in data_objects:
            data_file.append(f"/Positions/{particle}", buffers[particle])


def clear_buffer():
    buffers["Misc"] = pd.DataFrame(columns = ["time", "time-step", "distance between", "perturber_a", "inner_a"])
    for particle in data_objects:
        buffers[particle] = pd.DataFrame(columns = ["x", "y", "z"])


def save_data(sim):
    save_data_to_buffer(sim)
    if buffers["Misc"].shape[0] >= 10000:
        save_data_to_disk()
        clear_buffer()

def dump_record(record):
    with open("outcome.json", "w") as file:
        json.dump(record, file, indent=4)
    save_data_to_disk()

def save_final_data(outcome_record, sim, start_time):
    outcome_record["Time"] = time.time() - start_time
    outcome_record["Final Inner Eccentricity"] = sim.particles['inner'].calculate_orbit(primary=sim.particles['SMBH']).e
    outcome_record["Final Inner Semimaj"] = sim.particles['inner'].calculate_orbit(primary=sim.particles['SMBH']).a
    outcome_record["Final Perturber Eccentricity"] = sim.particles['perturber'].calculate_orbit(primary=sim.particles['SMBH']).e
    outcome_record["Final Perturber Semimaj"] = sim.particles['perturber'].calculate_orbit(primary=sim.particles['SMBH']).a
    outcome_record["Final SMBH xyz"] = sim.particles['SMBH'].xyz
    outcome_record["Final SMBH vxyz"] = sim.particles['SMBH'].vxyz
    outcome_record["Final Inner xyz"] = sim.particles['inner'].xyz
    outcome_record["Final Inner vxyz"] = sim.particles['inner'].vxyz
    outcome_record["Final Perturber xyz"] = sim.particles['perturber'].xyz
    outcome_record["Final Perturber vxyz"] = sim.particles['perturber'].vxyz

class UnboundException(Exception):
    pass


class CollisionException(Exception):
    pass


#####################################################################################

outcome_record = {"Time": 0, "Result": "Unfinished", "Formed Binary": False, # one of these outcome jsons exist for each simulation
                  "Minimum Separation (AU)": np.inf, "Minimum Separation (Rg)": np.inf, "Minimum Separation Step": 0, 
                  "Inner xyz at Min Sep": [0,0,0], "Inner vxyz at Min Sep": [0,0,0],
                  "Perturber xyz at Min Sep": [0,0,0], "Perturber vxyz at Min Sep": [0,0,0],
                  #"Minimum t_GW": np.inf, "Minimum t_GW Time": np.inf, 
                  #"Inner xyz at Min t_GW": [0,0,0], "Inner vxyz at Min t_GW": [0,0,0],
                  #"Perturber xyz at Min t_GW": [0,0,0], "Perturber vxyz at Min t_GW": [0,0,0],
                  "SMBH_m": 0, "inner_m": 0, "perturber_m" : 0,
                  "inner_M": 0, "perturber_M": 0, 
                  "inner_inc": 0, "perturber_inc": 0,
                  "inner_e": 0, "perturber_e": 0,
                  "inner_a": 0, "perturber_a": 0,
                  "Final Inner Eccentricity": 0, "Final Inner Semimaj": 0,
                  "Final Perturber Eccentricity": 0, "Final Perturber Semimaj": 0,
                  "Final SMBH xyz" : [0,0,0], 'Final SMBH vxyz' : [0,0,0],
                  "Final Inner xyz" : [0,0,0], 'Final Inner vxyz' : [0,0,0],
                  "Final Perturber xyz" : [0,0,0], 'Final Perturber vxyz' : [0,0,0],}

total_time_steps_completed = 0

buffers = {
    "Misc": pd.DataFrame(columns = ["time", "time-step", "distance between", "perturber_a", "inner_a"]),
}
for particle in data_objects:
    buffers[particle] = pd.DataFrame(columns = ["x", "y", "z"])
#####################################################################################

