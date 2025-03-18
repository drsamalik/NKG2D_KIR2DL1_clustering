import os, sys, shutil, functools, hashlib, glob
import numpy as np
from  time import sleep
import random
from rand_NKG_KIR_mixing import disjoint_mixing
from rand_sppark import rand_init_species
from run_spparks import setup_run
debug_flag = False
ensembles =  1000
time_pts =40
Nx = 91
Ny = 91
Nz = 1
dim = (Nx,Ny,Nz)
in_pxl = 0.04
nkg2d = 53
kir = 11686
hla = 4213
lck = 32064
vav = 18351
shp = 49871
random.seed(0)
current_path = os.getcwd()
Ulbps = [20,40,100,150,200,300,400]
for ulbp in Ulbps:
    inner_path = f'{current_path}/Ulbp_{ulbp}'
    os.makedirs(inner_path, exist_ok=True)
    print(inner_path)
    for ensemble in range(ensembles):
        disjoint_mixing(ensemble,dim,in_pxl,inner_path)
        rand_init_species(ensemble,nkg2d,kir,ulbp,hla,lck,vav,shp,dim,inner_path)
        setup_run(ensemble,time_pts,Nx,Ny,Nz,inner_path)
        os.system(f"sbatch {inner_path}/Run_main/Run_Sim/run_sim_"+str(ensemble)+".sh")
print('Done')