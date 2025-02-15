from Bio.PDB import PDBList
from Bio.PDB import parse_pdb_header
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import pandas as pd
import mdtraj as md
import os
from tqdm import tqdm
import numpy as np
import warnings
warnings.filterwarnings("ignore")
plist = PDBList()
parser = PDBParser()
io = PDBIO()
df = pd.read_csv('Top8000-best_hom95_pdb_chain_list.csv')
if os.path.exists('Top8000_dihe_seq.csv'):
    output_df = pd.read_csv('Top8000_dihe_seq.csv', index_col=0)
else:
    output_df = pd.DataFrame(columns=['PDB', 'Chain', 'Res 1', 'Res 2', 'Res 3', 'Phi', 'Psi'])
file_skip = open('skipped_pdb.txt', 'a')
dir = 'PDBs'

for n, row in tqdm(df.iterrows()):
    #Parse df
    pdb_code = str(row['pdb_id'])
    chain = row['chain']
    chain_int = ord(chain) - 65
    pdb_file_path = f'PDBs/{pdb_code}.pdb'
    if pdb_code not in output_df['PDB'].values or chain not in output_df[output_df['PDB'] == pdb_code]['Chain'].values:
        #Dowload PDB
        if not os.path.exists(pdb_file_path):
            plist.retrieve_pdb_file(pdb_code, pdir=dir, file_format='pdb')
            if not os.path.exists(f'PDBs/pdb{pdb_code}.ent'):
                continue
            os.rename(f'PDBs/pdb{pdb_code}.ent', pdb_file_path)

        #Save just the chain of interest
        if not os.path.exists(f'PDB_chain/{pdb_code}_{chain}.pdb'):
            structure = parser.get_structure(pdb_code, pdb_file_path)
            model = structure[0]
            chain_pick = model[chain]
            io.set_structure(chain_pick)
            io.save(f'PDB_chain/{pdb_code}_{chain}.pdb')
    
        #Get sequence
        seq_all = {record.id: record.seq for record in SeqIO.parse(pdb_file_path, 'pdb-seqres')}
        seq = list(seq_all[f'{pdb_code.upper()}:{chain}'])
        #Get PDB header
        header = parse_pdb_header(f'{dir}/{pdb_code}.pdb')
    
        #Determine residue number for first and last residue
        pdb_file = open(f'PDB_chain/{pdb_code}_{chain}.pdb', 'r').readlines()
        first_line = pdb_file[0].split(' ')
        while '' in first_line:
            first_line.remove('')
        if len(first_line[4]) == 1:
            first = int(first_line[5].strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
        else:
            first = int(first_line[4].strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))

        #Determine missing residues
        if header['has_missing_residues']:
            fail = False
            #Determine which residues are missing
            miss_res = header['missing_residues']
            miss_res_list = []
            for entry in miss_res:
                if entry['chain'] == chain:
                    miss_res_list.append(entry['ssseq'])

            #Determine residue number for last residue
            res_list = []
            for i, line in enumerate(pdb_file):
                line_s = line.split(' ')
                if line_s[0] != 'ATOM':
                    if len(res_list) ==0:
                        fail = True
                        break
                    last = int(res_list[-1])
                    break
                while '' in line_s:
                    line_s.remove('')
                if len(line_s[4]) == 1:
                    res_list.append(line_s[5].strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
                else:
                    res_list.append(line_s[4].strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
                  
            if fail:
                file_skip.write(f'{pdb_code} {chain}\n')
                continue

            #Determine if non-terminal residues are missing
            miss_res_NT = []
            for res in miss_res_list:
                if res>first and res<last:
                    miss_res_NT.append(int(res)-first)
                elif res<first:
                    del seq[0]
                else:
                    del seq[-1]
            miss_res_NT = np.array(miss_res_NT, dtype=int)

        #Load with mdtraj
        traj = md.load(f'PDB_chain/{pdb_code}_{chain}.pdb')
        traj = traj.atom_slice(traj.topology.select('protein'))

        #Compute dihedrals
        index_phi, phi = md.compute_phi(traj)
        index_psi, psi = md.compute_psi(traj)

        #Reformat dihedrals
        phi = phi[0,:] * (180/np.pi)
        psi = psi[0,:] * (180/np.pi)
        for a, angle in enumerate(phi):
            if angle < 0:
                phi[a] = angle + 360
        for a, angle in enumerate(psi):
            if angle < -120:
                psi[a] = angle + 360

        #Format df for saving
        seq_1, seq_2, seq_3, phi_str, psi_str, resnum = [],[],[],[],[],[]
        if header['has_missing_residues'] and len(miss_res_NT) != 0:
            if (len(seq) != len(phi) + len(miss_res_NT) + 1) or (len(seq) != len(psi) + len(miss_res_NT) + 2):
                file_skip.write(f'{pdb_code} {chain}\n')
                continue
            offset = 1
            for i in range(1, len(seq)-1):
                if i in miss_res_NT:
                    offset += 1
                elif i in miss_res_NT - 1 or i in miss_res_NT + 1:
                    continue
                else:
                    seq_1.append(seq[i-1])
                    seq_2.append(seq[i])
                    seq_3.append(seq[i+1])
                    phi_str.append(phi[i-offset-1])
                    psi_str.append(psi[i+1-offset-1])
                    resnum.append(int(i+1))
        else:
            if (len(seq) != len(phi) + 1) or (len(phi) != len(psi)):
                file_skip.write(f'{pdb_code} {chain}\n')
                continue
            for i in range(1, len(phi)-1):
                seq_1.append(seq[i-1])
                seq_2.append(seq[i])
                seq_3.append(seq[i+1])
                phi_str.append(phi[i-1])
                psi_str.append(psi[i])
                resnum.append(int(i+1))

        temp_df = pd.DataFrame({'PDB': pdb_code, 'Chain': chain, 'Resnum': resnum, 'Res 1': seq_1, 'Res 2': seq_2, 'Res 3': seq_3, 'Phi': phi_str, 'Psi': psi_str})
        output_df = pd.concat([output_df, temp_df], ignore_index=True)
        output_df.to_csv('Top8000_dihe_seq.csv')
