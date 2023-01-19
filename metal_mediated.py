"""General fucntion of this file is to find out the interactions which are
mediated by metals"""

import numpy as np

import pandas as pd

from read_pdb_regex import read_pdb
from ressidues import METAL_MEDIATED, METALS, METAL_DISTANCES

def get_distance(line1, line3):
    """This function gets the distance between 2 pdb ATOM lines
    param line1, line3: list of a pdb line with one part of line as one element. Example ['ATOM','C','OP1'......]
    Return dist: float distance thus calculated"""

    x1 = float(line1[6])
    y1 = float(line1[7])
    z1 = float(line1[8])
    x2 = float(line3[6])
    y2 = float(line3[7])
    z2 = float(line3[8])
    a = np.array((x1, y1, z1))
    b = np.array((x2, y2, z2))
    dist = np.linalg.norm(a-b)

    return dist

def get_focussed_lines(pdb_file, chains):
    """This function in general provides with lines from the pdb file which 
    include either a metal or a acceptor atom, acceptor atoms are determined by METAL_MEDIATED 
    dictionary in ressidues file.
    param pdb_file: name or path of the pdb file
    param chains: list of chains. Example ['A','B']
    Return metal_atoms_lines, donor_atom_lines: lists of such lines"""

    parsed, lines = read_pdb(pdb_file)
    metal_atoms_lines = []
    if len(chains) == 1:
        chain1 = chains[0]
        chain2 = chain1
    else:
        chain1 = chains[0]
        chain2 = chains[1]
    donor_atom_lines = {}
    for chain in chains:
        donor_atom_lines[chain] = []
    
    for line in parsed:
        # print(line[2]=='O')
        if line[3] in METAL_MEDIATED.keys() and line[2] in METAL_MEDIATED[line[3]] and line[4] in chains:
            chain = line[4]
            donor_atom_lines[chain].append(line)
        if line[0] == 'HETATM' and line[2] in METALS:
            metal_atoms_lines.append(line)

    return metal_atoms_lines, donor_atom_lines

def make_pd_dataframe(final_interactions):
    """this function just converts final_interactions to a padatframe"""
    
    df_dict = {'Focus Atom':[], 'Focus Atom Res.Name':[], 'Focus Atom Res.Number':[], 'Focus Atom Chain':[],'Res. Name1':[], 'Res. Number 1':[], 'Chain 1':[],'Res. Name2':[], 'Res. Number 2':[], 'Chain 2':[],'Distance (Å)':[]}

    for instance in final_interactions.keys():
        df_dict['Focus Atom'].append(instance.split()[2])
        df_dict['Focus Atom Res.Name'].append(instance.split()[3])
        df_dict['Focus Atom Res.Number'].append(instance.split()[5])
        df_dict['Focus Atom Chain'].append(instance.split()[4])
        df_dict['Res. Name1'].append(final_interactions[instance][0][0][3])
        df_dict['Res. Number 1'].append(final_interactions[instance][0][0][5])
        df_dict['Chain 1'].append(final_interactions[instance][0][0][4])
        df_dict['Res. Name2'].append(final_interactions[instance][0][1][3])
        df_dict['Res. Number 2'].append(final_interactions[instance][0][1][5])
        df_dict['Chain 2'].append(final_interactions[instance][0][1][4])
        df_dict["Distance (Å)"].append(f"{'%.2f' % final_interactions[instance][1][1][0]}, {'%.2f' %final_interactions[instance][1][1][1]}")
    
    df = pd.DataFrame(df_dict)

    return df

def check_intercation_metal(pdb_file, chains):
    """This fucntion checks the interactions feasibility using METAL_DISTANCES 
    dictionary which can be found in the ressidues file.
    Param pdb_file and chains are of same structure as above.
    Return Final_interactions: dict of final thus gotten interactions. 
    Structure of final_interactions:

    [metal_atom_line:
                     [donorline1, donorline2],
                     [2,
                     [distance_of metal_line and donorline1, distance_of metal_line and donorline2]
                     ]
                     ]"""


    metal_atoms_lines, donor_atom_lines = get_focussed_lines(pdb_file, chains)
    interactions = {}
    for metal_line in metal_atoms_lines:
        
        metal_line1 = ' '.join(str(e) for e in metal_line)
        for chain in donor_atom_lines.keys():
            for donor_line in donor_atom_lines[chain]:
                dist = get_distance(metal_line, donor_line)
                if dist <= METAL_DISTANCES[metal_line[2]]:
                    # print(metal_line)
                    # print(donor_line)
                    if metal_line1 not in interactions.keys():
                        interactions[metal_line1] =([donor_line],[1,[dist]])
                    elif metal_line1 in interactions.keys() and interactions[metal_line1][0][0][4] != donor_line[4]:
                        interactions[metal_line1][0].append(donor_line)
                        interactions[metal_line1][1][0] = 2
                        interactions[metal_line1][1][1].append(dist)
                    elif metal_line1 in interactions.keys() and interactions[metal_line1][0][0][4] == donor_line[4] and interactions[metal_line1][0] != donor_line and interactions[metal_line1][0][0][5] != donor_line[5]:
                        interactions[metal_line1][0].append(donor_line)
                        interactions[metal_line1][1][0] = 2
                        interactions[metal_line1][1][1].append(dist)
                    elif metal_line1 in interactions.keys() and interactions[metal_line1][0][0][4] == donor_line[4] and interactions[metal_line1][0] != donor_line and interactions[metal_line1][0][0][5] != donor_line[5] and interactions[metal_line1][1][0] == 2:
                        old_metal_line = metal_line1
                        metal_line1 = f'{metal_line1}_'
                        metal_line2 = f'{metal_line1}_'
                        # print(metal_line1)
                        # print(donor_line)
                        interactions[metal_line1] =([donor_line,interactions[old_metal_line][0][0]],[2,[dist,interactions[old_metal_line][1][1][0]]])
                        interactions[metal_line2] =([donor_line,interactions[old_metal_line][0][1]],[2,[dist,interactions[old_metal_line][1][1][1]]])
                        # print(interactions[metal_line1])
                    print(interactions[metal_line1])
                    print()

    final_interactions = {}

    for key1 in interactions.keys():
        if interactions[key1][1][0] == 2:
            final_interactions[key1] = interactions[key1]

    # print(len(final_interactions.keys()))
    # print(final_interactions)
    if len(final_interactions.keys())>0:
        return make_pd_dataframe(final_interactions)
    else: 
        return None
    
# final = check_intercation_metal('1L2X.pdb',['A','A'] )

# water_o_atoms_lines, donor_atom_lines = get_focussed_lines('3u43.pdb', ['A','B'])
# final = check_intercation(water_o_atoms_lines, donor_atom_lines )

# for key in final.keys():
#     print(key)
#     print(final[key][0])
#     print(final[key][1])
#     print()
# final.to_csv('metal.csv')

