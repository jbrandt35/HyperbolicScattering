import os
import shutil
import numpy as np
from time import sleep
import sys
import datetime
import json

print('start')

migration_type = "drag"
migration_constant = 10 ** 5.5
gr_type = "1PN+2.5PN"
randomize_M = True
randomize_inc = True
randomize_ecc = True
randomize_mass = True
smbh_mass = 1e6
inner_mass, perturber_mass = 60, 10 # if random_mass = False
inner_mass_lower_bound, inner_mass_upper_bound = 10, 150 # if random_mass = True
perturber_mass_lower_bound, perturber_mass_upper_bound = 10, 150
inc_std = 0.01 # radians

perturber_separation_range = (5.0, 5.0, 0.25)
inner_semimaj_range = (1000, 1000, 10)
num_in_each = 10000
#num_in_each = 5

perturber_separation_range = np.array(list(range(int(100 * perturber_separation_range[0]), int(100 * (perturber_separation_range[1] + perturber_separation_range[2])), int(100 * perturber_separation_range[2])))) / 100
inner_semimaj_range = np.array(list(range(int(100 * inner_semimaj_range[0]), int(100 * (inner_semimaj_range[1] + inner_semimaj_range[2])), int(100 * inner_semimaj_range[2])))) / 100


def return_bash(script, path = "utils"):
    working_dir = os.getcwd()
    os.chdir(path)
    output_stream = os.popen("./" + script)
    output = output_stream.read().strip()
    output_stream.close()
    os.chdir(working_dir)
    return int(output)


def run_bash(script, options, path = "utils"):
    working_dir = os.getcwd()
    os.chdir(path)
    os.system(f"sh {script} {options}")
    os.chdir(working_dir)


def build_queue(inner_semimajs, perturber_separations):
    q = []
    for inner_semimaj in inner_semimajs:
        for perturber_separation in perturber_separations:
            for i in range(1, num_in_each + 1):
                label = f"{np.round(inner_semimaj,2)} {np.round(perturber_separation,2)} {i} {migration_type} {migration_constant} {gr_type} {randomize_M} {randomize_inc} {randomize_ecc} {randomize_mass} {smbh_mass} {inner_mass} {perturber_mass} {inner_mass_lower_bound} {inner_mass_upper_bound} {perturber_mass_lower_bound} {perturber_mass_upper_bound} {inc_std}"
                q.append(label)
    return q

def build_queue_alt(inner_semimajs, perturber_separation):
    q = []
    for inner_semimaj in inner_semimajs:
        for i in range(1, num_in_each + 1):
            label = f"{np.round(inner_semimaj,2)} {np.round(perturber_separation,2)} {i}"
            q.append(label)
    return q


def print_to_stdout(text):
    sys.stdout.write("\n" + text)

#'''
def write_info_file(perturber_separation_range, inner_semimaj_range, migration_type="None", migration_constant=0, gr_type="None", randomize_M = True, randomize_inc = True, randomize_ecc = True, randomize_mass = False):
    label_lines = ["Perturber Range", "Inner Range", "Migration Type", "Migration Constant", "Random Mass", "SMBH Mass", "Inner Mass", "Perturber Mass", "Inner Mass Range", "Perturber Mass Range", "GR Type", "Random Mean Anomaly", "Random Inclination", "Inclination Standard Deviation", "Random Eccentricity"]
    value_lines = [str(perturber_separation_range), str(inner_semimaj_range), str(migration_type), str(migration_constant), str(randomize_mass), str(smbh_mass), str(inner_mass), str(perturber_mass), str([inner_mass_lower_bound, inner_mass_upper_bound]), str([perturber_mass_lower_bound, perturber_mass_upper_bound]), str(gr_type), str(randomize_M), str(randomize_inc), str(inc_std), str(randomize_ecc)]

    with open('/storage/home/hhive1/jbrandt35/data/HyperbolicSims/runs/current_sims/00_info/setup.txt', 'w+') as f:
        for i,x in enumerate(label_lines):
            f.write(x)
            f.write(': ')
            f.write(value_lines[i])
            f.write('\n')
#'''

'''
def write_info_file(perturber_separation_range, inner_semimaj_range, migration_type="None", migration_constant=0, gr_type="None", randomize_M = True, randmize_inc = True, randomize_ecc = True, randomize_mass = False):
    info_dict = {"Perturber Range": perturber_separation_range, "Inner Range": inner_semimaj_range, "Migration Type": migration_type, "Migration Constant": migration_constant, "GR Type": gr_type, "Random Mean Anomaly": randomize_M, "Random Inclination": randomize_inc, "Random Eccentricity": randomize_ecc, "Random Mass Ratio": randomize_mass}
    json_object = json.dumps(info_dict, indent=4)
 
    with open('/storage/home/hhive1/jbrandt35/data/HyperbolicSims/runs/current_sims/00_info/setup.json', 'w+') as outfile:
        outfile.write(json_object)
#'''
'''
write_info_file(perturber_separation_range, inner_semimaj_range, migration_type=migration_type, migration_constant=migration_constant, gr_type=gr_type, randomize_M = randomize_M, randomize_inc = randomize_inc, randomize_ecc = randomize_ecc, randomize_mass = randomize_mass)
#'''

label_lines = ["Perturber Range", "Inner Range", "Migration Type", "Migration Constant", "Random Mass", "SMBH Mass", "Inner Mass", "Perturber Mass", "Inner Mass Range", "Perturber Mass Range", "GR Type", "Random Mean Anomaly", "Random Inclination", "Inclination Standard Deviation", "Random Eccentricity"]
value_lines = [str(perturber_separation_range), str(inner_semimaj_range), str(migration_type), str(migration_constant), str(randomize_mass), str(smbh_mass), str(inner_mass), str(perturber_mass), str([inner_mass_lower_bound, inner_mass_upper_bound]), str([perturber_mass_lower_bound, perturber_mass_upper_bound]), str(gr_type), str(randomize_M), str(randomize_inc), str(inc_std), str(randomize_ecc)]

with open('/storage/home/hhive1/jbrandt35/data/HyperbolicSims/runs/current_sims/00_info/setup.txt', 'w+') as f:
    for i,x in enumerate(label_lines):
        f.write(x)
        f.write(': ')
        f.write(value_lines[i])
        f.write('\n')

queue = build_queue(inner_semimaj_range, perturber_separation_range)
#queue = build_queue_alt(inner_semimaj_range, 5.0)

while len(queue) > 0:
    print_to_stdout(f"There are {len(queue)} jobs in the queue")
    number_of_jobs_submittable = 495 - return_bash("running_jobs.sh")
    print_to_stdout(f"There are {number_of_jobs_submittable} slots to submit \n")
    for i in range(0, number_of_jobs_submittable):
        try:
            job = queue.pop()
            run_bash("build.sh", job)
        except IndexError:
            print_to_stdout("The queue is finished!")
            break
    sleep(60)








