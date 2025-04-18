import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from find_cluster import get_clusters
import copy
def ploting(matrix,title):
    fig,ax = plt.subplots(1,1,figsize=(10,4))
    ax.imshow(matrix, cmap='gray', vmin=matrix.min(), vmax=matrix.max())
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{title}',fontsize=12)
    plt.show()
def Receptor_site_init(recptor,r_sites,clusters):
    for i in range(recptor):
        cluster = random.choice(clusters)
        r_sites[cluster[0],cluster[1]] = r_sites[cluster[0],cluster[1]]+1
    return r_sites
def NKG_sites(num_mol, mix_cluters):
    init_site = np.zeros_like(mix_cluters)
    nkg_bi_intensity = copy.deepcopy(mix_cluters)
    nkg_bi_intensity[nkg_bi_intensity==2]=0
    nkg_clusters = get_clusters(nkg_bi_intensity)
    nkg_clusters = sum(nkg_clusters,[])
    nkg_site = Receptor_site_init(num_mol,init_site,nkg_clusters)
    #ploting(nkg_site+nkg_bi_intensity,title=f'nkg_r_sites_{num_mol}')
    return nkg_site,nkg_bi_intensity
def KIR_sites(num_mol, mix_cluters):
    init_site = np.zeros_like(mix_cluters)
    kir_bi_intensity = copy.deepcopy(mix_cluters)
    kir_bi_intensity[kir_bi_intensity==1]=0
    kir_bi_intensity[kir_bi_intensity==2]=1
    kir_clusters = get_clusters(kir_bi_intensity)
    kir_clusters = sum(kir_clusters,[])
    kir_site = Receptor_site_init(num_mol,init_site,kir_clusters)
    #ploting(kir_site+kir_bi_intensity,title=f'kir_r_sites_{num_mol}')
    return kir_site,kir_bi_intensity
def Molecule_site(num_mol,mix_cluters,mol_type):
    site = np.zeros_like(mix_cluters)
    for i in range(num_mol):
        cluster = (random.randint(0, mix_cluters.shape[0]-1), random.randint(0, mix_cluters.shape[1]-1))
        site[cluster[0],cluster[1]] = site[cluster[0],cluster[1]]+1
    #ploting(site,title=f'{mol_type}_{num_mol}')
    return site


def rand_init_species(ensemble,nkg2d,kir,ulbp3,hla,lck,vav,shp,dim,current_path):
    Nx,Ny,Nz=dim[0],dim[1],dim[2]
    mix_cluters = np.loadtxt(f'{current_path}/Init_files/NKG_KIR_disjoint_mixing_iter_{ensemble}.txt').astype(np.int64)
    f2=open(f'{current_path}/Init_files/input_sppark.nkg2d_kir_ulbp3_hla.disjoint_{ensemble}',"w")
    f2.write('Site file written by dump sites 1 command\n\n')
    f2.write('3 dimension\n')
    f2.write(f'{Nx*Ny} sites\n')
    f2.write('id i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 values\n')
    f2.write(f'0 {Nx} xlo xhi\n')
    f2.write(f'0 {Ny} ylo yhi\n')
    f2.write(f'0 {Nz} zlo zhi\n\n')
    f2.write('Values\n\n')
    nkg_r_sites,nkg_bi_intensity = NKG_sites(nkg2d,mix_cluters)
    kir_r_sites,kir_bi_intensity = KIR_sites(kir,mix_cluters)

    nkg_l_sites = Molecule_site(ulbp3,mix_cluters,mol_type='nkg_l_sites')
    kir_l_sites = Molecule_site(hla,mix_cluters,mol_type='kir_l_sites')
    lck_sites = Molecule_site(lck,mix_cluters,mol_type='lck_sites')
    vav_sites = Molecule_site(vav,mix_cluters,mol_type='vav_sites')
    shp_sites = Molecule_site(shp,mix_cluters,mol_type='shp_sites')

    for j in range(mix_cluters.shape[0]):
        for i in range(mix_cluters.shape[1]):
            site_id = j*mix_cluters.shape[1]+i
            f2.write("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n"%(site_id+1,nkg_l_sites[j,i],0,i+1,j+1,0,nkg_r_sites[j,i],0,0,lck_sites[j,i],vav_sites[j,i],0,0,0,0,0,nkg_bi_intensity[j,i],kir_bi_intensity[j,i],kir_l_sites[j,i],kir_r_sites[j,i],0,0,shp_sites[j,i]))
    f2.close()
            
    

        






    
