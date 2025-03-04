import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib.pyplot as plt
import argparse
import warnings
import random
from itertools import combinations_with_replacement, product

def find_cluster_centroid(df):
    torsion_dist_array = np.zeros((len(df.index), len(df.index)))
    n=0
    for i, row_i in df.iterrows():
        m=0
        for j, row_j in df.iterrows():
            torsion_dist_array[n][m] = torsion_distance(row_i, row_j)
            m += 1
        n += 1
    tor_sum = np.zeros(len(df.index))
    for t in range(len(df.index)):
        tor_sum[t] = sum(torsion_dist_array[t][:])
    min_index = np.argsort(tor_sum)[0]
    cid = list(df['Identifier'].values)[min_index]
    return cid

def get_torsion(input_df):
    for pdb in tqdm(input_df['PDB'].unique()):
        if pdb not in tor_4mer_df['PDB'].unique():
            pdb_df = input_df[input_df['PDB'] == pdb]
            for chain in pdb_df['Chain'].unique():
                chain_df = pdb_df[pdb_df['Chain'] == chain]
                num_res = len(chain_df.index)
                for res in range(1, num_res-3):
                    sect_df = chain_df[chain_df['Resnum'].isin([res, res+1, res+2, res+3])]
                    seq = ''.join(sect_df['Res 1'].values)
                    simp_seq_list = []
                    for aa in seq:
                        if aa in ['I', 'T', 'V']:
                            simp_seq_list.append('B')
                        elif aa == 'P' or aa == 'G':
                            simp_seq_list.append(aa)
                        else:
                            simp_seq_list.append('L')
                    simp_seq = ''.join(simp_seq_list)
                    df = pd.DataFrame({'PDB': pdb, 'Chain': chain, 'Res Range': f'{res}-{res+3}', 'Sequence': seq, 'Simplified Sequence': simp_seq, 'Phi 1': sect_df.iloc[0, 6], 'Psi 1': sect_df.iloc[0, 7], 
                                   'Phi 2': sect_df.iloc[1, 6], 'Psi 2': sect_df.iloc[1, 7],
                                   'Phi 3': sect_df.iloc[2, 6], 'Psi 3': sect_df.iloc[2, 7],
                                   'Phi 4': sect_df.iloc[3, 6], 'Psi 4': sect_df.iloc[3, 7]}, index=[len(tor_4mer_df.index)])
                    tor_4mer_df = pd.concat([tor_4mer_df, df])
            tor_4mer_df.to_csv('4mer_torsion_raw.csv')

def angle_diff(angle1, angle2):
    diff = np.abs(angle1 - angle2) % 360
    return min(diff, 360 - diff)

def torsion_distance(row1, row2):
    sq_sum = 0
    for i in range(1, 5):
        sq_sum += (angle_diff(row1[f'Phi {i}'], row2[f'Phi {i}']))**2 + (angle_diff(row1[f'Psi {i}'], row2[f'Psi {i}']))**2
    return np.sqrt(sq_sum)    

parser = argparse.ArgumentParser(description = 'Cluster 4-mers Based on Torsional Distance and Sequence')
parser.add_argument('-n', type=int, required=True, help='Total # of Conformers Desired')
parser.add_argument('-c', type=int, required=True, help='# of Clusters for Each Selected Sequence')

args = parser.parse_args()
num_tot_clusters = args.n
num_cluster = args.c

top8000_dihe = pd.read_csv('../Top8000/Top8000_dihe_seq.csv', index_col=0)
if not os.path.exists('4mer_torsion_raw.csv'):
    tor_4mer_df = pd.DataFrame()
else:
    tor_4mer_df = pd.read_csv('4mer_torsion_raw.csv', index_col=0)

if len(tor_4mer_df['PDB'].unique()) != len(top8000_dihe['PDB'].unique()):
    for pdb in tqdm(top8000_dihe['PDB'].unique()):
        pdb_df = top8000_dihe[top8000_dihe['PDB'] == pdb]
        for chain in pdb_df['Chain'].unique():
            if len(tor_4mer_df.index) == 0 or len(tor_4mer_df[(tor_4mer_df['PDB'] == pdb) & (tor_4mer_df['Chain'] == chain)].index) == 0:
                res_range, seq, simp_seq, phi1, phi2, phi3, phi4, psi1, psi2, psi3, psi4 = [],[],[],[],[],[],[],[],[],[],[]
                chain_df = pdb_df[pdb_df['Chain'] == chain]
                residues = chain_df['Resnum'].to_list()
                for r in range(len(residues)-1):
                    res_range_str = residues[r:r+4]
                    res1 = int(res_range_str[0])
                    res4 = int(res_range_str[-1])
                    if res4 == res1+3:
                        df_4mer = chain_df[chain_df['Resnum'].isin(res_range_str)]
                        seq_r = ''.join([df_4mer.iloc[0,3], df_4mer.iloc[0,4], df_4mer.iloc[3,2], df_4mer.iloc[3,3]])
                        seq.append(seq_r)
                        simp_seq_list = []
                        for s in seq_r:
                            if s == 'P':
                                simp_seq_list.append('P')
                            elif s == 'G':
                                simp_seq_list.append('G')
                            elif s in ['I', 'T', 'V']:
                                simp_seq_list.append('B')
                            else:
                                simp_seq_list.append('L')
                        simp_seq.append(''.join(simp_seq_list))
                        phi1.append(df_4mer.iloc[0,5])
                        psi1.append(df_4mer.iloc[0,6])
                        phi2.append(df_4mer.iloc[1,5])
                        psi2.append(df_4mer.iloc[1,6])
                        phi3.append(df_4mer.iloc[2,5])
                        psi3.append(df_4mer.iloc[2,6])
                        phi4.append(df_4mer.iloc[3,5])
                        psi4.append(df_4mer.iloc[3,6])
                        res_range.append(f'{res1}-{res4}')
                df = pd.DataFrame({'PDB': pdb, 'Chain': chain, 'Res Range': res_range, 'Sequence': seq, 'Simplified Sequence': simp_seq, 'Phi 1': phi1, 'Psi 1': psi1, 'Phi 2': phi2, 'Psi 2': psi2, 'Phi 3': phi3, 'Psi 3': psi3, 'Phi 4': phi4, 'Psi 4': psi4})
                tor_4mer_df = pd.concat([tor_4mer_df, df])
                tor_4mer_df.to_csv('4mer_torsion_raw.csv')
else:
    print('Completed 4mer_torsion_raw.csv file input')

if os.path.exists('4mer_simp_seq_freq.csv'):
    simp_seq_freq = pd.read_csv('4mer_simp_seq_freq.csv', index_col=0)
else:
    unique_seq = [''.join(combination) for combination in product(tuple('LGPB'), repeat=4)]
    num_4mer = np.zeros(len(unique_seq))
    for s, seq in enumerate(unique_seq):
        num_4mer[s] = len(tor_4mer_df[tor_4mer_df['Simplified Sequence'] == seq].index)
    simp_seq_freq = pd.DataFrame({'Sequence': unique_seq, 'Frequency': num_4mer})
    simp_seq_freq.to_csv('4mer_simp_seq_freq.csv')

total_freq = sum(simp_seq_freq['Frequency'].values)

if os.path.exists(f'4mer_seq_num_{num_tot_clusters}.csv'):
    num_seq_add = pd.read_csv(f'4mer_seq_num_{num_tot_clusters}.csv', index_col=0)
else:
    sequence_list = list(simp_seq_freq[simp_seq_freq['Frequency'] > 0]['Sequence'].values)
    num_seq_all = np.zeros(len(sequence_list))
    for s, simp_seq in enumerate(sequence_list):
        freq = simp_seq_freq[simp_seq_freq['Sequence'] == simp_seq]['Frequency'].values[0]
        num_seq = int((num_tot_clusters/num_cluster) * (freq / total_freq))
        num_seq_all[s] = num_seq
    prev_index = np.nan
    while sum(num_seq_all) > num_tot_clusters/num_cluster:
        sorted_indices = np.argsort(num_seq_all)
        if prev_index != sorted_indices[-1]:
            prev_index = sorted_indices[-1]
        else:
            prev_index = sorted_indices[-2]
        num_seq_all[prev_index] -= 1
    
    #Add to the simplified sequence with the highest frequency of all 0 sequences
    freq_array = simp_seq_freq['Frequency'].values
    max_indices = np.argsort(freq_array)[::-1] #indices from highest to lowest
    while sum(num_seq_all) < num_tot_clusters/num_cluster:
        zero_indices = np.where(num_seq_all == 0)[0]
        for i in max_indices:
            if i in zero_indices:
                index = i
                break
        num_seq_all[index] += 1
    num_seq_add = pd.DataFrame({'Simplified Sequence': sequence_list, 'Num Seq': num_seq_all})
    num_seq_add.to_csv(f'4mer_seq_num_{num_tot_clusters}.csv')

if os.path.exists('4mer_torsion_cluster.csv') and os.path.exists('selected_4mers.csv'):
    cluster_4mer_df = pd.read_csv('4mer_torsion_cluster.csv', index_col=0)
    selected_conformers = pd.read_csv('selected_4mers.csv', index_col=0)
else:
    cluster_4mer_df = pd.DataFrame()
    selected_conformers = pd.DataFrame()

if not os.path.exists('HC'):
    os.mkdir('HC')

list_all_simp_seq = num_seq_add[num_seq_add['Num Seq'] > 0]['Simplified Sequence'].values
for simp_seq in tqdm(list_all_simp_seq):
    if len(selected_conformers) == 0 or simp_seq not in selected_conformers['Simplified Sequence'].values:
        simp_seq_df = tor_4mer_df[tor_4mer_df['Simplified Sequence'] == simp_seq] 
        num_seq = num_seq_add[num_seq_add['Simplified Sequence'] == simp_seq]['Num Seq'].values[0]
        seq_list = list(simp_seq_df['Sequence'].unique())
        if len(seq_list) < 2 or num_seq > len(seq_list):
            seq_select = []
            for seq in seq_list:
                seq_df = tor_4mer_df[tor_4mer_df['Sequence'] == seq]
                if len(seq_df.index) > num_cluster:
                    seq_select.append(seq)
        else:
            seq_list_freq = np.zeros(len(seq_list))
            for s, seq_test in enumerate(seq_list):
                seq_list_freq[s] = len(tor_4mer_df[tor_4mer_df['Sequence'] == seq_test])
            sorted_index_seq_list_freq = np.flip(np.argsort(seq_list_freq))
            seq_select = []
            for idx in sorted_index_seq_list_freq:
                if seq_list_freq[idx] > num_cluster:
                    seq_select.append(seq_list[idx])
                if len(seq_select) == num_seq:
                    break
        for seq in seq_select:
            seq_df = tor_4mer_df[tor_4mer_df['Sequence'] == seq]
            if len(seq_df.index) == 0:
                print(seq)
                exit()
            torsion_dist_array = np.zeros((len(seq_df.index), len(seq_df.index)))
            n=0
            for i, row_i in seq_df.iterrows():
                m=0
                for j, row_j in seq_df.iterrows():
                    torsion_dist_array[n][m] = torsion_distance(row_i, row_j)
                    m += 1
                n += 1
            linkage_matrix = linkage(torsion_dist_array, method='single') 

            plt.figure(figsize=(10, 5))
            dendrogram(linkage_matrix)
            plt.xlabel("Data Points")
            plt.ylabel("Distance")
            plt.title("Hierarchical Clustering Dendrogram")
            plt.savefig(f'HC/hierarch_cluster_{seq}.png')
            plt.close()
    
            clusters = fcluster(linkage_matrix, t=num_cluster, criterion='maxclust')
            seq_df['Cluster'] = clusters
            cluster_4mer_df = pd.concat([cluster_4mer_df, seq_df])
            seq_df['Identifier'] = np.arange(len(seq_df.index))

            unique_cluster_id = list(set(clusters))
            for clust in unique_cluster_id:
                cluster_df = seq_df[seq_df['Cluster'] == clust]
                cid = find_cluster_centroid(cluster_df)
                selected_conformers = pd.concat([selected_conformers, cluster_df[cluster_df['Identifier'] == cid]])
        cluster_4mer_df.reset_index(inplace=True, drop=True)
        selected_conformers.reset_index(inplace=True, drop=True)
        cluster_4mer_df.to_csv('4mer_torsion_cluster.csv')
        selected_conformers.to_csv('selected_4mers.csv')
    
