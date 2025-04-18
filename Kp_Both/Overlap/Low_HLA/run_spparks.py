import os
import numpy as np
def setup_run(x,time_pts, Nx, Ny, Nz,inner_path):
    num_CPUs = 1
    sim_path = f'{inner_path}/Run_main/Run_Sim'
    os.makedirs(sim_path,exist_ok=True)
    slurm_path = f'{inner_path}/Run_main/Slurm_Out'
    os.makedirs(slurm_path,exist_ok=True)
    out_path = f'{inner_path}/Out_files'
    os.makedirs(out_path,exist_ok=True)
    with open(f'{sim_path}/run_sim_{x}.sh','w') as fRun:
        fRun.write("#!/bin/sh\n")
        #fRun.write("#SBATCH --mail-user=Saeed.Ahmad@nationwidechildrens.org\n")
        fRun.write("#SBATCH --mail-type=FAIL,REQUEUE\n")
        fRun.write("#SBATCH --job-name=20dis_inb\n")
        fRun.write("#SBATCH --time=13:00:00\n") # 7-0:00:00
        fRun.write("##SBATCH --partition=himem\n") # 7-0:00:00
        fRun.write("#SBATCH --cpus-per-task="+str(num_CPUs)+"\n")
        fRun.write("#SBATCH --ntasks=1\n")
        fRun.write(f"#SBATCH -o {slurm_path}/slurm.%j.out # STDOUT\n")
        fRun.write(f"#SBATCH -e {slurm_path}/slurm.%j.err # STDERR\n")
        fRun.write("\n")
        fRun.write("module load GCC/7.3.0-2.30 OpenMPI/3.1.1\n")
        fInName = f"{sim_path}/in."+str(x)
        fRun.write("mpirun -np 1 ./spk_redsky < "+fInName+"\n")    
        with open(fInName,'w') as fIn:
            fIn.write("# SPPARKS \n")
            fIn.write("\n")
            fIn.write(f"log      {out_path}/log.spparks_"+str(x)+"\n")
            fIn.write("seed     "+str(np.random.randint(0, 10000)+1)+"\n")
            fIn.write("\n")
            fIn.write("app_style    inb/diff/custom\n")
            fIn.write("\n")
            fIn.write("dimension    3\n")
            fIn.write("lattice      sc/6n 1\n")
            fIn.write("region       box block 0 "+str(Nx)+" 0 "+str(Ny)+" 0 "+str(Nz)+"\n")
            fIn.write("create_box   box\n")
            fIn.write("create_sites box\n")
            fIn.write("\n")
            fIn.write("region       global_bottom block 0 91 0 91 0 0\n")
            fIn.write("\n")
            old_path = out_path.replace('NKR_no_kp', 'NKR_General')
            fIn.write(f"read_sites       {inner_path}/Init_files/input_sppark.nkg2d_kir_ulbp3_hla.overlap_"+str(x)+"\n")
            fIn.write("\n")
            fIn.write("set      i2 value 1 region global_bottom\n")
            fIn.write("\n")
            fIn.write("sector   yes\n")
            fIn.write("solve_style  tree \n")
            fIn.write("\n")
            fIn.write("# species (note that in the C++ code indexing starts at 0 instead of 1)\n")
            fIn.write("add_species  l # i1 - NKG2DL\n")
            fIn.write("add_species  E # i2 - empty space (i.e. not a wall or free surface)\n")
            fIn.write("add_species  y # i3 - X cord \n")
            fIn.write("add_species  Y # i4 - Y cord \n")
            fIn.write("add_species  N # i5 - Voxel No.a\n")
            fIn.write("add_species  a # i6 - NKG2D\n")
            fIn.write("add_species  A # i7 - NKG2D-NKG2DL\n")
            fIn.write("add_species  p # i8 - pNKG2D-NKG2DL\n")
            fIn.write("add_species  k # i9 - SFK kinases\n")
            fIn.write("add_species  v # i10 -  unbounb Vav1s\n")
            fIn.write("add_species  V # 11 -  bounb Vav1s\n")
            fIn.write("add_species  P # i12 - phosphorylated Vav1 attached to AR_p\n")
            fIn.write("add_species  c # i13 - Lck+A complex \n")
            fIn.write("add_species  C # i14 - AR_p +vav+bounb SFK complex \n")
            fIn.write("add_species  f # i15 - free/unbound phosphorylated vav\n")
            fIn.write("add_species  W # i16 - AR cluster\n")
            fIn.write("add_species  w # i17 - IR cluster\n")
            fIn.write("add_species  L # i18 - unbound inhibiting ligands\n")
            fIn.write("add_species  i # i19 - unbound inhibiting receprors\n")
            fIn.write("add_species  I # i20 - bound inhibiting receprors\n")
            fIn.write("add_species  h # i21 - phosphorylated bound inhibiting receprors\n")
            fIn.write("add_species  s # i22 - unbound SHP\n")
            fIn.write("add_species  S # i23 - bound SHP\n")
            fIn.write("add_species  x # i24 - Lck+I complex\n")
            fIn.write("add_species  G # i25 - SHP attached to IR_p + p_Vav attached to AR_p\n")
            fIn.write("add_species  g # i26 - SHP attached to IR_p + free p_Vav\n")
            fIn.write("add_species  b # i27 - global sim box\n")
            fIn.write("add_species  B # i28 - boundary of sim_box (next to wall)\n")
            fIn.write("add_species  r # i29 - unbound vav'\n")
            fIn.write("add_species  R # i30 - bound vav'\n")
            fIn.write("add_species  D # i31 - ARp+bound vav'\n")
            fIn.write("add_species  Q # i32 - phos vav'+ARp\n")
            fIn.write("add_species  F # i33 - free phos vav'\n")
            fIn.write("add_species  Z # i34 - Shp+phos vav'+ARp\n")
            fIn.write("add_species  z # i35 - SHP+free phos vav'\n")
            fIn.write("\n")
            fIn.write("# NKG2D + NKG2DL1 <--> NKG2D-NKG2DL\n")
            fIn.write("add_rxn      0 local l a nbr 12.3854688 local A nbr\n")
            fIn.write("add_rxn      1 local A nbr 0.023 local l a nbr\n")
            fIn.write("\n")
            fIn.write("# diffusion NKG2DL, NKG2D and SFKs\n")
            fIn.write("\n")
            fIn.write("add_rxn      2 local l nbr E 6.25 local nbr l E\n")
            fIn.write("add_rxn      3 local a nbr W 6.25 local nbr a W\n")
            fIn.write("add_rxn      4 local k nbr E 6.25 local nbr k E\n")
            fIn.write("\n")
            fIn.write("# NKG2D-NKG2DL + SFK <--> NKG2D-NKG2DL-SFK\n")
            fIn.write("add_rxn      5 local A k nbr 50.0 local c nbr\n")
            fIn.write("add_rxn      6 local c nbr 0.006 local A k nbr\n")
            fIn.write("\n")
            fIn.write("# NKG2D-NKG2DL-SFK --> pNKG2D-NKG2DL + SFK\n")
            fIn.write("add_rxn      7 local c nbr 2.1899058353952916 local p k nbr\n")
            fIn.write("\n")
            fIn.write("# pNKG2D-NKG2DL + Vav <--> pNKG2D-NKG2DL-Vav\n")
            fIn.write("add_rxn      8 local p v nbr 0.63421875 local V nbr\n")
            fIn.write("add_rxn      9 local V nbr 0.01 local p v nbr\n")
            fIn.write("\n")
            fIn.write("# pNKG2D-NKG2DL-Vav + SFK <--> pNKG2D-NKG2DL-Vav-SFK\n")
            fIn.write("add_rxn      10 local V k nbr 50.0 local C nbr\n")
            fIn.write("add_rxn      11 local C nbr 0.028317078876353527 local V k nbr\n")
            fIn.write("\n")
            fIn.write("# pNKG2D-NKG2DL-Vav-SFK --> pNKG2D-NKG2DL-pVav +SFK\n")
            fIn.write("add_rxn      12 local C nbr 0.7762549321707577 local P k nbr\n")
            fIn.write("\n")
            fIn.write("# dephospho\n")
            fIn.write("add_rxn      13 local P nbr 1.0 local V nbr\n")
            fIn.write("add_rxn      14 local p nbr 2.0 local A nbr\n")
            fIn.write("\n")
            fIn.write("# dissociation complex-phosphoVav to AR_p + phosp-vav\n")
            fIn.write("add_rxn      15 local P nbr 0.01 local p f nbr\n")
            fIn.write("add_rxn      16 local p f nbr 0.63421875 local P nbr\n")
            fIn.write("\n")
            fIn.write("# dephospho of free pvav\n")
            fIn.write("add_rxn      17 local f nbr 1.0 local v nbr\n")
            fIn.write("\n")
            fIn.write("# dissociation of Ar complexes \n")
            fIn.write("add_rxn      18 local c nbr 0.023 local l a k nbr\n")
            fIn.write("add_rxn      19 local p nbr 0.023 local l a nbr\n")
            fIn.write("add_rxn      20 local V nbr 0.023 local l a v nbr\n")
            fIn.write("add_rxn      21 local C nbr 0.023 local l a v k nbr\n")
            fIn.write("add_rxn      22 local P nbr 0.023 local l a v nbr\n")
            fIn.write("\n")
            fIn.write("# diffusion\n")
            fIn.write("add_rxn      23 local i nbr w 6.25 local nbr i w\n")
            fIn.write("\n")
            fIn.write("# IR + L_I\n")
            fIn.write("add_rxn      24 local L i nbr 50.0 local I nbr\n")
            fIn.write("add_rxn      25 local I nbr 1.0 local L i nbr\n")
            fIn.write("\n")
            fIn.write("# AR+kinase\n")
            fIn.write("add_rxn      26 local I k nbr 50.0 local x nbr\n")
            fIn.write("add_rxn      27 local x nbr 0.006 local I k nbr\n")
            fIn.write("\n")
            fIn.write("# phophorylation of  IR\n")
            fIn.write("add_rxn      28 local x nbr 2.1899058353952916 local h k nbr\n")
            fIn.write("\n")
            fIn.write("# IR_p + Shp\n")
            fIn.write("add_rxn      29 local h s nbr 0.375 local S nbr\n")
            fIn.write("add_rxn      30 local S nbr 0.000509 local h s nbr\n")
            fIn.write("\n")
            fIn.write("# dephospho of IRp\n")
            fIn.write("add_rxn      31 local h nbr 2.0 local I nbr\n")
            fIn.write("\n")
            fIn.write("# (bound_SHP + IR_p complex) + (p_VAv + AR_p) complex\n")
            fIn.write("add_rxn      32 local S P nbr 50.0 local G nbr\n")
            fIn.write("add_rxn      33 local G nbr 0.01 local S P nbr\n")
            fIn.write("\n")
            fIn.write("# dephospho of dephosphorylation of pVAVcomplex due to bound-SHP\n")
            fIn.write("add_rxn      34 local G nbr 10.0 local S V nbr\n")
            fIn.write("\n")
            fIn.write("# (bound_SHP + IR_p complex) + (p_VAv + AR_p) complex\n")
            fIn.write("add_rxn      35 local S f nbr 6.3421875 local g nbr\n")
            fIn.write("add_rxn      36 local g nbr 0.01 local S f nbr\n")
            fIn.write("\n")
            fIn.write("# dephospho of dephosphorylation of pVAVcomplex due to bound-SHP\n")
            fIn.write("add_rxn      37 local g nbr 10.0 local S v nbr\n")
            fIn.write("\n")
            fIn.write("# dissociation of IR complexes \n")
            fIn.write("add_rxn      38 local x nbr 1.0 local L i k nbr\n")
            fIn.write("add_rxn      39 local h nbr 1.0 local L i nbr\n")
            fIn.write("add_rxn      40 local S nbr 1.0 local L i s nbr\n")
            fIn.write("add_rxn      41 local G nbr 1.0 local L i s P nbr\n")
            fIn.write("add_rxn      42 local G nbr 0.023 local l a v S nbr\n")
            fIn.write("add_rxn      43 local g nbr 1.0 local L i s v nbr\n")
            fIn.write("\n")
            fIn.write("# diffusion\n")
            fIn.write("add_rxn      44 local L nbr E 6.25 local nbr L E\n")
            fIn.write("\n")
            fIn.write("\n")
            fIn.write("diag_style   propensity\n")
            fIn.write("diag_style   array i1 sum i2 sum i3 sum i4 sum i5 sum i6 sum i7 sum i8 sum i9 sum i10 sum i11 sum i12 sum i13 sum i14 sum i15 sum i16 sum i17 sum i18 sum i19 sum i20 sum i21 sum i22 sum i23 sum i24 sum i25 sum i26 sum i27 sum i28 sum i29 sum i30 sum i31 sum i32 sum i33 sum i34 sum i35 sum\n")
            fIn.write("\n")
            fIn.write("stats           1.0\n")
            fIn.write(f"dump            1 sites 1.0 {out_path}/sites."+str(x)+".* id i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 i23 i24 i25 i26 i27 i28 i29 i30 i31 i32 i33 i34 i35\n") 
            fIn.write("\n")
            fIn.write("run             "+str(time_pts)+"\n")
