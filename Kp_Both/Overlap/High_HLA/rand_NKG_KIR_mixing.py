import numpy as np
import matplotlib.pyplot as plt
import copy,os
from matplotlib.colors import ListedColormap
from find_cluster import NKG_Cluster, KIR_Cluster, get_clusters
from com_relocation import random_relocation
from pathlib import Path
def disjoint_mixing(ensemble,dim,in_pxl,current_path):
    input_path = Path(current_path).parents[2]/'NKG_KIR_Init_data'
    NKG_Clusters, NKG_new_size, NKG_bi_intensity = NKG_Cluster(dim,in_pxl,input_path)
    KIR_Clusters, KIR_new_size, KIR_bi_intensity = KIR_Cluster(dim,in_pxl,input_path)
    NKG_relocated_cluster = random_relocation(NKG_Clusters,NKG_new_size,dim[:-1],nucleated_matrix=None)
    total_clusters=(NKG_new_size+KIR_new_size)
    Mix_cluters = random_relocation(KIR_Clusters,total_clusters,dim[:-1],nucleated_matrix=copy.deepcopy(NKG_relocated_cluster))
    os.makedirs(f'{current_path}/Init_files', exist_ok=True)
    np.savetxt(f'{current_path}/Init_files/NKG_KIR_overlap_mixing_iter_{ensemble}.txt', Mix_cluters, fmt='%d')
    #export_sites(f'NKG_KIR_disjoint_mixing_iter_{ensemble}.{Mix_cluters.shape[1]}.{Mix_cluters.shape[0]}', Mix_cluters)
    #cluster_plot(NKG_bi_intensity,KIR_bi_intensity,Mix_cluters,ensemble)

def export_sites(filename, binary_matrix):
    with open(filename, 'w') as f:
        # Write header
        f.write("Site file written by dump sites 1 command\n\n")
        f.write("3 dimension\n")
        f.write(f"{binary_matrix.shape[0]*binary_matrix.shape[1]} sites\n")
        f.write("id i3 i4 i6 i7 values\n")
        f.write(f"0 {binary_matrix.shape[1]} xlo xhi\n")
        f.write(f"0 {binary_matrix.shape[0]} ylo yhi\n")
        f.write(f"0 1 zlo zhi\n\n")
        f.write("Values\n\n")
        count = 0
        for j in range(binary_matrix.shape[0]):
            for i in range(binary_matrix.shape[1]):
                count += 1
                if int(binary_matrix[j,i])==0:
                    NKG,KIR = 0,0
                elif int(binary_matrix[j,i])==1:
                    NKG,KIR = 1,0
                else:
                    NKG,KIR = 0,1
                f.write(f"{count}\t{i+1}\t{j+1}\t{NKG}\t{KIR}\n")

def cluster_plot(NKG_bi_intensity,KIR_bi_intensity,Mix_cluters,ensemble):
    cmap0 = ListedColormap(['#000000', '#00FF00'])
    cmap1 = ListedColormap(['#000000', '#FF0000'])
    cmap2 = ListedColormap(['#000000', '#00FF00', '#FF0000'])
    fig,ax = plt.subplots(1,3,figsize=(21,6))
    ax[0].imshow(NKG_bi_intensity, cmap=cmap0, vmin=0, vmax=1)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title('NKG Intensity',fontsize=20)
    ax[1].imshow(KIR_bi_intensity, cmap=cmap1, vmin=0, vmax=1)
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_title('KIR Intensity',fontsize=20)
    ax[2].imshow(Mix_cluters, cmap=cmap2, vmin=0, vmax=2)
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    ax[2].set_title('Relocated Mix Intensity',fontsize=20)
    cbar = plt.colorbar(ax[2].get_images()[0], ax=ax[2], ticks=[0.35, 1, 1.70], shrink=0.95)
    cbar.ax.set_yticklabels(['0', 'NKG', 'KIR'],fontsize=20)
    plt.tight_layout()