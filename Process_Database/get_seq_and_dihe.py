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

def get_sequence_from_pdb(pdb_file):
    """
    Extracts the protein sequence from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        seq (str): The protein sequence as a string
    """
    res_dict = {'ARG': 'R', 'HIS': 'H', 'LYS': 'K', 'ASP': 'D', 'GLU': 'E', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'CYS': 'C', 'SEC': 'U', 'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W'}
    seq = ''
    res_num = None
    res_num_list = []
    for line in open(pdb_file, 'r').readlines():
        if line[:4] == 'ATOM':
            if line[17:20] == 'HOH':
                break
            if res_num is None or res_num != line[22:27]:
                if line[17:20] not in res_dict.keys():
                    res_code = 'X'
                else:
                    res_code = res_dict[line[17:20]]
                seq += res_code
                res_num = line[22:27]
                res_num_list.append(line[22:26])
    res_num_list_reset = np.zeros(len(res_num_list))
    modifier = 1
    for r, res in enumerate(res_num_list):
        if res == res_num_list[r-1]:
            modifier += 1
        res_num_list_reset[r] = int(res) - int(res_num_list[0]) + modifier
    return seq, res_num_list_reset

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
                file_skip.write(f'{pdb_code} {chain} File Not Found\n')
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
        full_seq = list(seq_all[f'{pdb_code.upper()}:{chain}'])
        seq, res_num_list = get_sequence_from_pdb(f'PDB_chain/{pdb_code}_{chain}.pdb')

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
            #Determine which residues are missing
            miss_res = header['missing_residues']
            miss_res_list = []
            for entry in miss_res:
                if entry['chain'] == chain:
                    miss_res_list.append(entry['ssseq'])

        else:
            res_num_list = np.arange(1, len(seq)+1, step=1)

        #Load with mdtraj
        traj = md.load(f'PDB_chain/{pdb_code}_{chain}.pdb')
        traj = traj.atom_slice(traj.topology.select(f'protein and resid < {len(res_num_list)}'))

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
        if (len(seq) != len(phi) + 1):
#            traj.save('test.pdb')
#            print(traj)
#            print(f'{pdb_code} {chain}\n')
#            print(len(seq))
#            print(seq)
#            print(len(phi))
#            exit()
            file_skip.write(f'{pdb_code} {chain}\n')
            continue
        
        for i in range(1, len(psi)-1):
            if res_num_list[i] == res_num_list[i-1]+1 and res_num_list[i] == res_num_list[i+1]-1:
                seq_1.append(seq[i-1])
                seq_2.append(seq[i])
                seq_3.append(seq[i+1])
                phi_str.append(phi[i-1])
                psi_str.append(psi[i])
                resnum.append(res_num_list[i])

        temp_df = pd.DataFrame({'PDB': pdb_code, 'Chain': chain, 'Resnum': resnum, 'Res 1': seq_1, 'Res 2': seq_2, 'Res 3': seq_3, 'Phi': phi_str, 'Psi': psi_str})
        output_df = pd.concat([output_df, temp_df], ignore_index=True)
        output_df.to_csv('Top8000_dihe_seq.csv')
