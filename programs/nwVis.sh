#!/bin/bash
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=50Gb
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --job-name=edgeList_creation

echo "STARTING at $(date)"

cd /home/pgsb/jan.ramos/programs

python3 pyfile.py

echo "ENDING at $(date)"

