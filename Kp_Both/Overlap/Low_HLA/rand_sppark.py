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
def ploting_clusters(all_clusters,occuped_cluters,title):
    fig,ax = plt.subplots(1,2,figsize=(10,4))
    ax[0].imshow(all_clusters, cmap='gray', vmin=all_clusters.min(), vmax=all_clusters.max())
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title(f'{title[0]}',fontsize=12)
    ax[1].imshow(occuped_cluters, cmap='gray', vmin=occuped_cluters.min(), vmax=occuped_cluters.max())
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title(f'{title[1]}',fontsize=12)
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
    nkg_clusters0 = get_clusters(nkg_bi_intensity)
    nkg_clusters1 = sum(nkg_clusters0,[])
    nkg_site = Receptor_site_init(num_mol,init_site,nkg_clusters1)
    occupied_nkg_cluster = copy.deepcopy(init_site)
    nkg_clusters0 = [cluster for cluster in nkg_clusters0 if any(nkg_site[y, x] > 0 for y, x in cluster)]
    occupied_nkg_cluster[tuple(np.array(sum(nkg_clusters0, [])).T)] = 1
    #ploting_clusters(nkg_site+nkg_bi_intensity,nkg_site+occupied_nkg_cluster,title=[f'{num_mol}_nkg_r_sites_in_{a}',f'{num_mol}_nkg_occupied_sites_in_{len(nkg_clusters0)}'])
    return nkg_site,occupied_nkg_cluster
def KIR_sites(num_mol, mix_cluters,nkg_bi_intensity):
    init_site = np.zeros_like(mix_cluters)
    occupied_kir_cluster  = np.zeros_like(mix_cluters)
    kir_bi_intensity = copy.deepcopy(mix_cluters)
    kir_bi_intensity[kir_bi_intensity==1]=0
    kir_bi_intensity[kir_bi_intensity==2]=1
    kir_overlapping_clusters = get_clusters(nkg_bi_intensity)
    kir_total_clusters = get_clusters(kir_bi_intensity)
    kir_overlapping_molecules = int((num_mol/len(kir_total_clusters))*len(kir_overlapping_clusters))
    kir_overlapping_sites = Receptor_site_init(kir_overlapping_molecules,init_site,sum(copy.deepcopy(kir_overlapping_clusters),[]))
    kir_nonovelapping_cluster = random.sample(kir_total_clusters,int(len(kir_total_clusters)-len(kir_overlapping_clusters)))
    kir_all_sites = Receptor_site_init(num_mol-kir_overlapping_molecules,copy.deepcopy(kir_overlapping_sites),sum(copy.deepcopy(kir_nonovelapping_cluster),[]))
    kir_all_cluster = kir_nonovelapping_cluster+kir_overlapping_clusters
    occupied_kir_cluster[tuple(np.array(sum(kir_all_cluster, [])).T)] = 1
    #ploting_clusters(kir_all_sites+occupied_kir_cluster,kir_overlapping_sites+nkg_bi_intensity,title=[f'kir_r_sites_{kir_all_sites.sum()}',f'kir_overlapping_sites_{kir_overlapping_sites.sum()}'])
    return kir_all_sites,occupied_kir_cluster
def Molecule_site(num_mol,mix_cluters,mol_type):
    site = np.zeros_like(mix_cluters)
    for i in range(num_mol):
        cluster = (random.randint(0, mix_cluters.shape[0]-1), random.randint(0, mix_cluters.shape[1]-1))
        site[cluster[0],cluster[1]] = site[cluster[0],cluster[1]]+1
    #ploting(site,title=f'{mol_type}_{num_mol}')
    return site


def rand_init_species(ensemble,nkg2d,kir,ulbp3,hla,lck,vav,shp,dim,current_path):
    Nx,Ny,Nz=dim[0],dim[1],dim[2]
    mix_cluters = np.loadtxt(f'{current_path}/Init_files/NKG_KIR_overlap_mixing_iter_{ensemble}.txt').astype(np.int64)
    f2=open(f'{current_path}/Init_files/input_sppark.nkg2d_kir_ulbp3_hla.overlap_{ensemble}',"w")
    f2.write('Site file written by dump sites 1 command\n\n')
    f2.write('3 dimension\n')
    f2.write(f'{Nx*Ny} sites\n')
    f2.write('id i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20 i21 i22 values\n')
    f2.write(f'0 {Nx} xlo xhi\n')
    f2.write(f'0 {Ny} ylo yhi\n')
    f2.write(f'0 {Nz} zlo zhi\n\n')
    f2.write('Values\n\n')
    nkg_r_sites,nkg_bi_intensity = NKG_sites(nkg2d,mix_cluters)
    kir_r_sites,kir_bi_intensity = KIR_sites(kir,mix_cluters,nkg_bi_intensity)
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
            
    

        






    
