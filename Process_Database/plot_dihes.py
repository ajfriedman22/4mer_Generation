import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_scatterplot(data, file_name, hue_cat=None, palette=None, add_guide=False):
    fig, ax = plt.subplots()
    if hue_cat != None:
        sns.scatterplot(data = data, x='Phi', y='Psi', hue=hue_cat, s=0.8, alpha=0.5)
    else:
        sns.scatterplot(data = data, x='Phi', y='Psi', s=0.8)
    if add_guide:
        ax.add_patch(Rectangle((180,105), 90, 90, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\beta$", (220,140), fontsize=18)
        ax.add_patch(Rectangle((60,-90), 45, 60, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\gamma$", (75,-65), fontsize=18)
        ax.add_patch(Rectangle((225,-60), 90, 105, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\delta$", (260,-10), fontsize=18)
        ax.add_patch(Rectangle((285,-60), 30, 30, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\alpha$", (290,-50), fontsize=18)
        ax.add_patch(Rectangle((45,120), 135, 120, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\varepsilon$", (105,180), fontsize=18)
        ax.add_patch(Rectangle((195,45), 60, 60, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\zeta$", (220,70), fontsize=18)
        ax.add_patch(Rectangle((270,120), 45, 75, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$P_{II}$", (285,145), fontsize=18)
        ax.add_patch(Rectangle((255,45), 45, 60, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\gamma'$", (270,70), fontsize=18)
        ax.add_patch(Rectangle((30,-15), 90, 90, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$\delta'$", (70,25), fontsize=18)
        ax.add_patch(Rectangle((45,165), 75, 75, edgecolor='black', fill=False, lw=3))
        plt.annotate(r"$P_{II}'$", (80,200), fontsize=18)
    plt.xlabel('Phi(deg)')
    plt.ylabel('Psi(deg)')
    plt.xlim(0,360)
    plt.ylim(-120,240)
    if hue_cat != None:
        plt.legend(loc='upper right')
    fig.savefig(f'{file_name}.png')
    plt.close(fig)

def plot_heatmap(data, file_name, deg_int, add_guide=False, diff=False):
    fig, ax = plt.subplots()
    if diff:
        sns.heatmap(data, cmap="vlag", center=0, xticklabels=np.linspace(10, 350, num=int(360/deg_int), dtype=int), yticklabels=np.linspace(-110, 230, num=int(360/deg_int), dtype=int))
    else:
        sns.heatmap(data, cmap="YlOrBr", xticklabels=np.linspace(10, 350, num=int(360/deg_int), dtype=int), yticklabels=np.linspace(-110, 230, num=int(360/deg_int), dtype=int))
    if add_guide:
        ax.add_patch(Rectangle((int(180/deg_int),int((105+120)/deg_int)), 90/deg_int, 90/deg_int, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\beta$", (int(220/deg_int),int((140+120)/deg_int)), fontsize=16)
        ax.add_patch(Rectangle((int(225/deg_int),int((-60+120)/deg_int)), 90/deg_int, 105/deg_int, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\delta$", (int(260/deg_int),int((-10+120)/deg_int)), fontsize=16)
        ax.add_patch(Rectangle((int(285/deg_int),int((-60+120)/deg_int)), 30/deg_int, 30/deg_int, edgecolor='black', fill=False, lw=2))
        plt.annotate(r"$\alpha$", (int(290/deg_int),int((-50+120)/deg_int)), fontsize=18)
    ax.invert_yaxis()
    plt.xlabel('Phi Bin')
    plt.ylabel('Psi Bin')
    fig.savefig(file_name)
    plt.close(fig)

def reformat_data(df):
    for i, row in df.iterrows():
        if int(row['Phi']) < 0:
            df.at[i,'Phi'] = row['Phi'] + 360
        if int(row['Psi']) < -120:
            df.at[i, 'Psi'] = row['Psi'] + 360
    return df

if not os.path.exists('dihe_top8000_processed.csv'):
    top8000 = pd.read_csv('Top8000_dihe_seq.csv')
    print('Data Loaded')

    #Reformat data
    top8000_processed = reformat_data(top8000)
    top8000_processed.to_csv('dihe_top8000_processed.csv')
    print('Processed data')
else:
    top8000_processed = pd.read_csv('dihe_top8000_processed.csv')
    print('Processed Data Loaded')

#Plot dihedrals
plot_scatterplot(top8000_processed, 'Top8000_dihe_w_guide', add_guide=True)

#Create dihedral bins and determine the weight of those bins
deg_int=15
deg_bin_top8000 = np.zeros((int(360/deg_int), int(360/deg_int))) #Phi, Psi
phi_bin, psi_bin = [], []
for index, row in top8000_processed.iterrows():
    phi = row['Phi']
    psi = row['Psi']
    pdb = row['PDB']
    i = int(phi/deg_int)
    j = int((psi+120)/deg_int)
    deg_bin_top8000[j][i] += 1
    phi_bin.append(i)
    psi_bin.append(j)
top8000_df = top8000_processed[['PDB', 'Chain', 'Res 1', 'Res 2', 'Res 3']].copy()
top8000_df['Phi Bin'] = phi_bin
top8000_df['Psi Bin'] = psi_bin
    
deg_bin_top8000_per = deg_bin_top8000/len(top8000_processed.index)

plot_heatmap(deg_bin_top8000_per, 'Top8000_dihe_bin', 15, add_guide=True)

