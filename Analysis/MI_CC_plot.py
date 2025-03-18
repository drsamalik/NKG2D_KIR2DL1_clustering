import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
def all_plos(matrix,cluster_condn,color,plot_info=None):
    categories = ['kp','no kp inh','no kp act','no kp']
    subcategories = ['Low ','High']
    values = matrix
    colors = [color]*len(categories)
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['font.family'] = 'Arial'
    size = 12
    ax0 = plt.gca()
    ax0.tick_params(axis='both', which='major', labelsize=19, length=5, width=3)
    bar_width = 0.5
    x = np.arange(len(categories))*1.5  # Base positions for categories
    offset = np.array([-bar_width / 1.8, bar_width / 2])  # Offsets for subcategories
    #patterns = [['..', '.O'],['**','*O'],['|./','|o-'],['\\\o\||','--..||']]  # Different fill patterns for categories
    patterns = [['|./','|o-'],['/','/'],['\\','\\'],['O','O']]  # Different fill patterns for categories
    
    # Plotting the bars
    for i, (cat, vals, cols) in enumerate(zip(categories, values, colors)):
        for j, (val, col) in enumerate(zip(vals, cols)):
            if i ==0:
                ax.bar(x[i] + offset[j], val, width=bar_width-0.03, color=col, edgecolor='black',linewidth=2)  
            else:
                ax.bar(x[i] + offset[j], val, width=bar_width-0.03, color=col, edgecolor='black',linewidth=2, hatch=patterns[i][j])  
            ax.text(x[i] + offset[j], val+0.01, f'{val:.2f}', ha='center', va='bottom', fontsize=size, fontweight='bold')  
    
    
    
    # Formatting the axes
    # ax.set_ylabel(plot_info, fontsize=20, fontweight='bold')
    ax.set_xticks([-0.25,0.25,1.25,1.75,2.75,3.25,4.25,4.75])
    ax.set_xticklabels('', fontsize=size, fontweight='bold')
    ax.set_yticks([0.0,0.8,1.6])
    # Adding the subcategory labels under the x-axis
    # for i, cat in enumerate(categories):
    #     ax.text(x[i] - 1.0*bar_width / 2, -0.06, subcategories[0], ha='center', va='top', fontsize=size, fontweight='bold')
    #     ax.text(x[i] + 1.15*bar_width / 2, -0.06, subcategories[1], ha='center', va='top', fontsize=size, fontweight='bold')
    #     ax.text(x[i] - bar_width / 2, -0.15, 'HLA', ha='center', va='top', fontsize=size, fontweight='bold')
    #     ax.text(x[i] + 1.15*bar_width / 2, -0.15, 'HLA', ha='center', va='top', fontsize=size, fontweight='bold')
    #     #ax.text(x[i] , -0.21, '------------', ha='center', va='top', fontsize=size, fontweight='bold')
    #     #ax.text(x[i], -0.215, '[', ha='center', va='top', fontsize=30, fontweight='light', rotation=90)
    #     ax.text(x[i], -0.3, cat, ha='center', va='top', fontsize=size, fontweight='bold', style='italic')
    plt.ylim(0,1.76)
    # plt.title(cluster_condn, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{plot_info}_plot_{cluster_condn}.png')