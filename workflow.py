import os
import sys

import ressidues as RESIDUES
from inter_non_covalent_with_all_possibilities import make_files_for_all_rings
from segregation_all_possibilites import creating_sub_pdbs
from orient_non_covalent_refactor import read_sub_pdb
from filter_refactor import filtering
from Plotting import plotting

def magic_begin():
    pdb_file = sys.argv[1]
    chain_1 = sys.argv[2]
    chain_2 = sys.argv[3]
    threshold = sys.argv[4]
    rings = RESIDUES.all_ring_names
    output_file_name = f'{os.path.split(pdb_file)[1]}_feasible_ring_interactions.txt'
    output_file_with_interaction = open(output_file_name, 'w')
    chains = [chain_1,chain_2]

    new_dir = f'{os.path.split(pdb_file)[1]}_sub_pdb'
    os.mkdir(new_dir)
    os.mkdir(f'{new_dir}/Reoriented')
    os.mkdir(f'{new_dir}/Reoriented/xyz')
    os.mkdir(f'{new_dir}/Reoriented/data')
    os.mkdir(f'{new_dir}/Reoriented/filtered')
    print(output_file_name, pdb_file)

    make_files_for_all_rings(rings,pdb_file, chains, output_file_with_interaction)
    creating_sub_pdbs(output_file_name,pdb_file,new_dir)
    read_sub_pdb(new_dir)
    filtering(new_dir,'cation_py')
    plotting(pdb_file, chains, threshold,f'{new_dir}/Reoriented/filtered')
# reading_all_file(output_file_name)