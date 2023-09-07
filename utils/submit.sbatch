#!/bin/bash

#SBATCH -JHypNAME
#SBATCH --account=hive-gli323
#SBATCH -t24:00:00
#SBATCH -phive
#SBATCH -ooutput.out

cd WORKDIR
module load anaconda3/2021.05
conda activate /storage/home/hhive1/jbrandt35/.conda/envs/myenv
srun python /storage/home/hhive1/jbrandt35/data/HyperbolicSims/utils/rare_use/text_to_json.py
srun python HyperbolicSimulation.py
srun python HyperbolicPostProc.py
rm *.py