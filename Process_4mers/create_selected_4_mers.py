import mdtraj as md
import numpy as np
import pandas as pd
import os
from openbabel import openbabel
from tqdm import tqdm
from Bio import SeqIO
from Bio.PDB import PDBParser
import re
from Bio.PDB import PDBList
from Bio.PDB import parse_pdb_header
from Bio.PDB.PDBIO import PDBIO

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
    for line in open(pdb_file, 'r').readlines():
        if line[:4] == 'ATOM':
            if line[17:20] == 'HOH':
                break
            if res_num is None or res_num != line[22:26]:
                if line[17:20] not in res_dict.keys():
                    res_code = 'X'
                else:
                    res_code = res_dict[line[17:20]]
                seq += res_code
                res_num = line[22:26]

    return seq

def process_line(line):
    if len(''.join(filter(str.isdigit, line[4]))) != 0:
        line.insert(4, 'A')
        line[5] = ''.join(filter(str.isdigit, line[5]))
    if len(line[6]) > 8:
        text = line[6]
        decimal_loc = re.findall(r"\d+\.\d+", text)
        try:
            end = decimal_loc[0]+3
        except:
            print(line[6])
            print(decimal_loc)
        line[6] = text[:end]
        line.insert(7, text[end+1:])
    if len(line[7]) > 8:
        text = line[7]
        decimal_loc = text.index('.')
        end = decimal_loc+4
        line[7] = text[:end]
        line.insert(8, text[end:])
    return line

if not os.path.exists('selected_4mers'):
    os.mkdir('selected_4mers')
if not os.path.exists('selected_6res'):
    os.mkdir('selected_6res')
plist = PDBList()
parser = PDBParser()
io = PDBIO()
dir = 'PDBs'
selected_4mers = pd.read_csv('selected_4mers.csv', index_col=0)

#Add MDTraj formatted residue ranges
mdtraj_input = []
for i, row in selected_4mers.iterrows():
    pdb = row['PDB']
    chain = row['Chain']
    res1, res2 = row['Res Range'].split('-')
    pdb_file = open(f'../Top8000/PDB_chain/{pdb}_{chain}.pdb')
    first_res = None
    l=0
    prev_res = None
    for line in pdb_file:
        if (line[0:4] == 'ATOM' or line[0:4] == 'HETA') and line[13:15] == 'CA' and line[23:27] != prev_res:
            prev_res = line[23:27]
            if first_res is None:
                first_res = int(line[22:26])
                if first_res == -2:
                    first_res = 1
            if int(line[22:26]) - first_res + 1 == int(res1):
                start = l
                break
            l += 1
    try:
        mdtraj_input.append(f'resid >= {start-1} and resid <= {start+4}')
    except:
        print(f'{pdb} - {chain} - {first_res} - {res1}')
    del start
        
phi1, phi2, phi3, phi4, psi1, psi2, psi3, psi4 = [],[],[],[],[],[],[],[]
error_index = []
dihe_err = pd.DataFrame()
for seq in tqdm(selected_4mers['Sequence'].unique()):
    if os.path.exists(f'selected_4mers/{seq}.sdf'):
        continue
    if not os.path.exists(f'selected_6res/{seq}'):
        os.mkdir(f'selected_6res/{seq}')
    if not os.path.exists(f'selected_4mers/{seq}'):
        os.mkdir(f'selected_4mers/{seq}')
    seq_df = selected_4mers[selected_4mers['Sequence'] == seq]
    for i, row in seq_df.iterrows():
        #Create a 6-mer to use the extra residues to set the dihedrals
        res_range = mdtraj_input[i]
        pdb = row['PDB']
        chain = row['Chain']
        if not os.path.exists(f'selected_6res/{seq}/{i}.pdb'):
            full_struct = md.load(f'../Top8000/PDB_chain/{pdb}_{chain}.pdb')
            struct_6res = full_struct.atom_slice(full_struct.topology.select(res_range))
            if struct_6res.n_residues != 6:
                struct_6res.save('test.pdb')
                raise Exception(f'Error on row {i}')
            struct_6res.save(f'selected_6res/{seq}/{i}.pdb')
            #Check that sequence matches
            seq_from_4mer = get_sequence_from_pdb(f'selected_6res/{seq}/{i}.pdb')
            if seq_from_4mer[1:-1] != row['Sequence']:
                error_index.append(i)

        #Reformat the terminal residues to look like ACE and NME
        if not os.path.exists(f'selected_4mers/{seq}/{i}_noH.pdb'):
            init_pdb = open(f'selected_6res/{seq}/{i}.pdb', 'r').readlines()
            output_pdb = open(f'selected_4mers/{seq}/{i}_noH.pdb', 'w')
            output_pdb.write(init_pdb[0])
            atom_id, res_num = 1, 0
            prev_res = None
            for line in init_pdb[1:-2]:
                if line[:4] == 'ATOM':
                    if line[13:15] == 'N ':
                        res_num += 1
                    if '\n' in line:
                        end = ''
                    else:
                        end = '\n'
                    if line[77] == 'H':
                        continue
                    if res_num == 1: #Change to ACE
                        if line[13:15] not in ['C ', 'O ', 'CA']:
                            continue
                        output_pdb.write(line[:4] + str(atom_id).rjust(7) + '  ' + line[13:16].ljust(4) + 'ACE' + 'A'.rjust(2) + str(res_num).rjust(4) + line[26:] + end)
                    elif res_num == 6: #Change to NME
                        if line[13:15] not in ['N ', 'CA']:
                            continue
                        output_pdb.write(line[:4] + str(atom_id).rjust(7) + '  ' + line[13:16].ljust(4) + 'NME' + 'A'.rjust(2) + str(res_num).rjust(4) + line[26:] + end)
                    else: #Just write
                        output_pdb.write(line[:4] + str(atom_id).rjust(7) + '  ' + line[13:16].ljust(4) + line[17:20] + 'A'.rjust(2) + str(res_num).rjust(4) + line[26:] + end)
                    atom_id += 1
                else:
                    output_pdb.write(line[:4] + str(atom_id).rjust(8) + 'NME'.rjust(9) + 'A'.rjust(2) + str(res_num).rjust(4) + end)
            output_pdb.write(init_pdb[-2]) 
            output_pdb.write(init_pdb[-1])
            output_pdb.close()
    
        #Add hydrogens and save in sdf format
        if not os.path.exists(f'selected_4mers/{seq}/{i}.sdf'):
            # Initialize the Open Babel conversion object
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("pdb", "sdf")

            # Create an empty molecule
            mol = openbabel.OBMol()

            # Read the PDB file
            obConversion.ReadFile(mol, f'selected_4mers/{seq}/{i}_noH.pdb')

            # Add hydrogens
            mol.AddHydrogens()

            # Write the modified molecule to a new PDB file
            obConversion.WriteFile(mol, f'selected_4mers/{seq}/{i}.sdf')

        #Check Dihedrals
        traj = md.load(f'selected_4mers/{seq}/{i}_noH.pdb')
        index, phi = md.compute_phi(traj)
        index, psi = md.compute_psi(traj)
        phi = phi[0,:] * (180/np.pi)
        psi = psi[0,:] * (180/np.pi)
        for a, angle in enumerate(phi):
            if angle < 0:
                phi[a] = angle + 360
        for a, angle in enumerate(psi):
            if angle < -120:
                psi[a] = angle + 360
        add_dict = {'PDB': pdb, 'Chain': chain, 'Residue Range': res_range}
        for r in range(4):
            add_dict[f'Phi {r+1}'] = phi[r]
            add_dict[f'Psi {r+1}'] = psi[r]
        dihe_err = pd.concat([dihe_err, pd.DataFrame(add_dict, index=[len(dihe_err)])])
    #Combine sdf files
    if not os.path.exists(f'selected_4mers/{seq}.sdf'):
        atom_in_each_file = []
        output_file = open(f'selected_4mers/{seq}.sdf', 'w')
        for i in seq_df.index:
            add_file = open(f'selected_4mers/{seq}/{i}.sdf', 'r').readlines()
            atoms_in_file = []
            for line in add_file:
                output_file.write(line)
                if len(line) > 50:
                    atoms_in_file.append(line[31])
            atom_in_each_file.append(atoms_in_file)
        output_file.close()

        #Check that atoms match exactly
        if len(atom_in_each_file[0][:]) != len(atom_in_each_file[1][:]) or len(atom_in_each_file[0][:]) != len(atom_in_each_file[2][:]) or len(atom_in_each_file[0][:]) != len(atom_in_each_file[3][:]) or len(atom_in_each_file[0][:]) != len(atom_in_each_file[4][:]):
            raise Exception(f'Atoms not matching in {seq}')
        
        for a, atom in enumerate(atom_in_each_file[0]):
            if atom != atom_in_each_file[1][a] or atom != atom_in_each_file[2][a] or atom != atom_in_each_file[3][a] or atom != atom_in_each_file[4][a]:
                raise Exception(f'Atoms not matching in {seq}')

error_df = selected_4mers.iloc[error_index]
error_df.to_csv('error_seq.csv')
dihe_err.to_csv('final_4mer_dihe.csv')
