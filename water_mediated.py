import numpy as np

import pandas as pd
from read_pdb_regex import read_pdb
from ressidues import WATER_MEDIATED

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
    include either a water or a donor atom, donor atoms are determined by WATER_MEDIATED 
    dictionary in ressidues file.
    param pdb_file: name or path of the pdb file
    param chains: list of chains. Example ['A','B']
    Return metal_atoms_lines, donor_atom_lines: lists of such lines"""

    parsed, lines = read_pdb(pdb_file)
    water_o_atoms_lines = []
    if len(chains) == 1:
        chain1 = chains[0]
        chain2 = chain1
    else:
        chain1 = chains[0]
        chain2 = chains[1]
    donor_atom_lines = {}
    for chain in chains:
        donor_atom_lines[chain] = []
    # donor_atom_lines = {chain1:[], chain2:[]}
    for line in parsed:
        # print(line[2]=='O')
        if line[3] in WATER_MEDIATED.keys() and line[2] in WATER_MEDIATED[line[3]] and line[4] in chains:
            chain = line[4]

            # print(line)
            donor_atom_lines[chain].append(line)
        if line[0] == 'HETATM' and line[2] == 'O':
            water_o_atoms_lines.append(line)

    return water_o_atoms_lines, donor_atom_lines

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
        df_dict['Distance (Å)'].append(f"{'%.2f' % final_interactions[instance][1][1][0]}, {'%.2f' %final_interactions[instance][1][1][1]}")
    
    df = pd.DataFrame(df_dict)

    return df


def check_intercation_water(pdb_file, chains):
    """This fucntion checks the interactions feasibility using distance of 3.5 Å
    Param pdb_file and chains are of same structure as above.
    Return Final_interactions: dict of final thus gotten interactions. 
    Structure of final_interactions:

    [water_atom_line:
                     [donorline1, donorline2],
                     [2,
                     [distance_of water_line and donorline1, distance_of water_line and donorline2]
                     ]
                     ]"""

    water_o_atoms_lines, donor_atom_lines = get_focussed_lines(pdb_file, chains)
    interactions = {}
    for o_line in water_o_atoms_lines:
        o_line1 = ' '.join(str(e) for e in o_line)
        for chain in donor_atom_lines.keys():
            for donor_line in donor_atom_lines[chain]:
                dist = get_distance(o_line, donor_line)
                if dist <= 3.5:
                    # print(o_line)
                    # print(donor_line)
                    # print(dist)
                    if o_line1 not in interactions.keys():
                        interactions[o_line1] =([donor_line],[1,[dist]])
                    elif o_line1 in interactions.keys() and interactions[o_line1][0][0][4] != donor_line[4]:
                        interactions[o_line1][0].append(donor_line)
                        interactions[o_line1][1][0] = 2
                        interactions[o_line1][1][1].append(dist)
                    elif o_line1 in interactions.keys() and interactions[o_line1][0][0][4] == donor_line[4] and interactions[o_line1][0] != donor_line and interactions[o_line1][0][0][5] != donor_line[5]:
                        interactions[o_line1][0].append(donor_line)
                        interactions[o_line1][1][0] = 2
                        interactions[o_line1][1][1].append(dist)
                    
    final_interactions = {}

    for key1 in interactions.keys():
        if interactions[key1][1][0] == 2:
            final_interactions[key1] = interactions[key1]

    """ Final_interactions = {O atom line as string:
                            [[first chain interactor, second chain interactor],
                            [2,[distance_first_chain, distance_second_chain]]}"""
    if len(final_interactions.keys())>0:
        return make_pd_dataframe(final_interactions)
    else: 
        return None


# final = check_intercation_water('1L2X.pdb',['A'] )

# for key in final.keys():
#     print(key)
#     print(final[key][0][0])
#     print(final[key][0][1])
#     print(final[key][1][1])
#     print()

# print(final)
# final.to_csv('water.csv')


# O--OP1 (2.96)   o::OP2 (3.17)