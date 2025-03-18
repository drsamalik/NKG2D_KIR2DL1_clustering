#!/bin/sh
#SBATCH --mail-user=Saeed.Ahmad@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --job-name=20dis_inb
#SBATCH --time=13:00:00
##SBATCH --partition=himem
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH -o /gpfs0/home1/gddaslab/sxa126/NKR_Nano_Clustering/NKR_NKG2D_Codes/No_Kp/Overlap/Low_HLA/Ulbp_20/Run_main/Slurm_Out/slurm.%j.out # STDOUT
#SBATCH -e /gpfs0/home1/gddaslab/sxa126/NKR_Nano_Clustering/NKR_NKG2D_Codes/No_Kp/Overlap/Low_HLA/Ulbp_20/Run_main/Slurm_Out/slurm.%j.err # STDERR

module load GCC/7.3.0-2.30 OpenMPI/3.1.1
mpirun -np 1 ./spk_redsky < /gpfs0/home1/gddaslab/sxa126/NKR_Nano_Clustering/NKR_NKG2D_Codes/No_Kp/Overlap/Low_HLA/Ulbp_20/Run_main/Run_Sim/in.0
