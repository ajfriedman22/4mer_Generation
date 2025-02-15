import mdtraj as md
import numpy as np
import pandas as pd
import os
from openbabel import openbabel
from tqdm import tqdm
from Bio import SeqIO
from Bio.PDB import PDBParser
import re

def get_sequence_from_pdb(pdb_file):
    """
    Extracts the protein sequence from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        str: The protein sequence as a string, or None if an error occurs.
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        
        # Iterate through models and chains, assuming only one protein chain
        for model in structure:
            for chain in model:
                # Extract sequence from SEQRES records if available
                for record in SeqIO.parse(pdb_file, "pdb-seqres"):
                    if record.id.split(":")[1] == chain.id:
                        return str(record.seq)
                # If SEQRES is not found, extract from ATOM records (less reliable)
                sequence = ""
                for residue in chain:
                    sequence += residue.get_resname()[:1]  # Use first letter of residue name
                return sequence
    except Exception as e:
        print(f"Error processing PDB file: {e}")
        return None

def process_line(line):
    if len(''.join(filter(str.isdigit, line[4]))) != 0:
        line.insert(4, 'A')
        line[5] = ''.join(filter(str.isdigit, line[5]))
    if len(line[6]) > 8:
        text = line[6]
        decimal_loc = re.findall(r"\d+\.\d+", text)
        end = decimal_loc[0]+3
        line[6] = text[:end]
        line.insert(7, text[end+1:])
    if len(line[7]) > 8:
        text = line[7]
        decimal_loc = text.index('.')
        end = decimal_loc+4
        line[7] = text[:end]
        line.insert(8, text[end:])
    return line


selected_4mers = pd.read_csv('selected_4mers.csv', index_col=0)
error_index = []
phi1, phi2, phi3, phi4, psi1, psi2, psi3, psi4 = [],[],[],[],[],[],[],[]
for i, row in tqdm(selected_4mers.iterrows()):
    #Create a 6-mer to use the extra residues to set the dihedrals
    res_range = row['Res Range']
    res1, res2 = res_range.split('-')
    if not os.path.exists(f'selected_6res/{i}.pdb'):
        pdb = row['PDB']
        chain = row['Chain']
        full_struct = md.load(f'../Top8000/PDB_chain/{pdb}_{chain}.pdb')
        if int(res1) == 1:
            struct_6res = full_struct.atom_slice(full_struct.topology.select(f'resid < {int(res2)+2}'))
        elif int(res1) == 2:
            struct_6res = full_struct.atom_slice(full_struct.topology.select(f'resid < {int(res2)+1}'))
        else:
            struct_6res = full_struct.atom_slice(full_struct.topology.select(f'resid > {int(res1)-3} and resid < {int(res2)+1}'))
        if struct_6res.n_residues != 6:
            raise Exception(f'Error on row {i}')
        struct_6res.save(f'selected_6res/{i}.pdb')
    #Check that sequence matches
    seq = get_sequence_from_pdb(f'selected_6res/{i}.pdb')
    if seq[1:-1] != row['Sequence']:
        error_index.append(i)
        continue

    #Reformat the terminal residues to look like ACE and NME
    if not os.path.exists(f'selected_4mers/{i}_noH.pdb'):
        init_pdb = open(f'selected_6res/{i}.pdb', 'r').readlines()
        output_pdb = open(f'selected_4mers/{i}_noH.pdb', 'w')
        output_pdb.write(init_pdb[0])
        atom_id, res_num = 1, 0
        prev_res = None
        for line in init_pdb[1:-2]:
            split_line = line.split(' ')
            while '' in split_line:
                split_line.remove('')
            if split_line[0] == 'ATOM':
                if split_line[2] == 'N':
                    res_num += 1
                if len(split_line) != 13:
                    split_line = process_line(split_line)
                n = 6
                if '\n' in split_line[n+5]:
                    end = ''
                else:
                    end = '\n'
                if res_num == 1: #Change to ACE
                    if split_line[2] not in ['C', 'O', 'CA']:
                        continue
                    output_pdb.write(split_line[0] + str(atom_id).rjust(7) + '  ' + split_line[2].ljust(4) + 'ACE' + 'A'.rjust(2) + str(res_num).rjust(4) + split_line[n].rjust(12) + split_line[n+1].rjust(8) + split_line[n+2].rjust(8) + split_line[n+3].rjust(6) + split_line[n+4].rjust(10) + split_line[n+5].rjust(8) + end)
                elif res_num == 6: #Change to NME
                    if split_line[2] not in ['N', 'CA']:
                        continue
                    output_pdb.write(split_line[0] + str(atom_id).rjust(7) + '  ' + split_line[2].ljust(4) + 'NME' + 'A'.rjust(2) + str(res_num).rjust(4) + split_line[n].rjust(12) + split_line[n+1].rjust(8) + split_line[n+2].rjust(8) + split_line[n+3].rjust(6) + split_line[n+4].rjust(10) + split_line[n+5].rjust(8) + end)
                else: #Just write
                    output_pdb.write(split_line[0] + str(atom_id).rjust(7) + '  ' + split_line[2].ljust(4) + split_line[3] + 'A'.rjust(2) + str(res_num).rjust(4) + split_line[n].rjust(12) + split_line[n+1].rjust(8) + split_line[n+2].rjust(8) + split_line[n+3].rjust(6) + split_line[n+4].rjust(10) + split_line[n+5].rjust(8) + end)
                atom_id += 1
            else:
                output_pdb.write(split_line[0] + str(atom_id).rjust(8) + 'NME'.rjust(9) + 'A'.rjust(2) + str(res_num).rjust(4) + end)
        output_pdb.write(init_pdb[-2]) 
        output_pdb.write(init_pdb[-1])
        output_pdb.close()
    
    #Add hydrogens and save in sdf format
    if not os.path.exists(f'selected_4mers/{i}.sdf'):
        # Initialize the Open Babel conversion object
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "sdf")

        # Create an empty molecule
        mol = openbabel.OBMol()

        # Read the PDB file
        obConversion.ReadFile(mol, f'selected_4mers/{i}_noH.pdb')

        # Add hydrogens
        mol.AddHydrogens()

        # Write the modified molecule to a new PDB file
        obConversion.WriteFile(mol, f'selected_4mers/{i}.sdf')

    #Check Dihedrals
    #traj = md.load(f'selected_4mers/{i}_noH.pdb')
    #index, phi = md.compute_phi(traj)
    #index, psi = md.compute_psi(traj)
    #phi = phi[0,:] * (180/np.pi)
    #psi = psi[0,:] * (180/np.pi)
    #for a, angle in enumerate(phi):
    #    if angle < 0:
    #        phi[a] = angle + 360
    #for a, angle in enumerate(psi):
    #    if angle < -120:
    #        psi[a] = angle + 360

error_df = selected_4mers.iloc[error_index]
error_df.to_csv('error_seq.csv')

