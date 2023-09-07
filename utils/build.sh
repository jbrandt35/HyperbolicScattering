
inner_semimaj=$1
perturber_separation=$2
i=$3
migration_type=$4
migration_constant=$5
gr_type=$6
randomize_M=$7
randomize_inc=$8
randomize_ecc=$9
randomize_mass=${10}
smbh_mass=${11}
inner_mass=${12}
perturber_mass=${13}
inner_mass_lower_bound=${14}
inner_mass_upper_bound=${15}
perturber_mass_lower_bound=${16}
perturber_mass_upper_bound=${17}
inc_std=${18}

mkdir -p ../runs/current_sims/run_data/${i}
mkdir -p ../runs/current_sims/run_data/${i}/plots

cp code/*.py ../runs/current_sims/run_data/${i}
cp submit.sbatch ../runs/current_sims/run_data/${i}

cd ../runs

sed -i "s|NAME|${i}|g" current_sims/run_data/${i}/submit.sbatch
sed -i "s|WORKDIR|/storage/home/hhive1/jbrandt35/data/HyperbolicSims/runs/current_sims/run_data/${i}|g" current_sims/run_data/${i}/submit.sbatch

sed -i "s|PERTSEPARATION|${perturber_separation}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|INNERSEMIMAJ|${inner_semimaj}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|MIGRATION_TYPE|${migration_type}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|MIGRATION_CONST|${migration_constant}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|GR_TYPE|${gr_type}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|RANDOMIZE_MASS|${randomize_mass}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|RANDOMIZE_M|${randomize_M}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|RANDOMIZE_INC|${randomize_inc}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|INCLINATION_STD|${inc_std}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|RANDOMIZE_ECC|${randomize_ecc}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|SMBH_MASS|${smbh_mass}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|INNER_MASS_LOWER_BOUND|${inner_mass_lower_bound}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|INNER_MASS_UPPER_BOUND|${inner_mass_upper_bound}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|PERTURBER_MASS_LOWER_BOUND|${perturber_mass_lower_bound}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|PERTURBER_MASS_UPPER_BOUND|${perturber_mass_upper_bound}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|INNER_MASS|${inner_mass}|g" current_sims/run_data/${i}/HyperbolicSimulation.py
sed -i "s|PERTURBER_MASS|${perturber_mass}|g" current_sims/run_data/${i}/HyperbolicSimulation.py

sbatch current_sims/run_data/${i}/submit.sbatch

