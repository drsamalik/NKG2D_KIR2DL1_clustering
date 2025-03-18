#!/bin/sh
#SBATCH --mail-user=Rajdeep.Grewal@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --job-name=main_20dis
#SBATCH --time=23:00:00
##SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
##SBATCH --exclude=gpu[01-05]
#SBATCH -o /home/gddaslab/sxa126/NKR_Nano_Clustering/NKR_Previous_Tries/NKR_General/Homogeneous/KMC_High_HLA/slurm.%j.out # STDOUT
#SBATCH -e /home/gddaslab/sxa126/NKR_Nano_Clustering/NKR_Previous_Tries/NKR_General/Homogeneous/KMC_High_HLA/slurm.%j.err # STDERR

echo 'Loading modules'
module purge 
module load miniforge3/24.3.0
conda activate PT_2.3.1_CUDA_12.1
echo 'activated conda environment'
python -u run_ensemble.py
conda deactivate
echo 'Ending modules'