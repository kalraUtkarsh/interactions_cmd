import numpy as np
import math

import pandas as pd
import vg

from read_pdb_regex import read_pdb
from ressidues import HALOGEN_ELECTRONEGATIVE, RADII, HALOGENS


def get_distance(line1, line3):
    x1 = float(line1[6])
    y1 = float(line1[7])
    z1 = float(line1[8])
    x2 = float(line3[6])
    y2 = float(line3[7])
    z2 = float(line3[8])

    a = np.array((x1, y1, z1))
    b = np.array((x2, y2, z2))

    dist = np.linalg.norm(a - b)

    return dist


# def dvivide_in_residues(parsed, lines):


def get_focussed_lines(pdb_file, chains):
    parsed, lines = read_pdb(pdb_file)
    halogen_atoms_lines = []
    chain1 = chains[0]
    chain2 = chains[1]
    donor_atom_lines = {chain1: [], chain2: []}
    residue_chunks = {}

    for line in parsed:
        # print(line[2]=='O')
        if line[5] not in residue_chunks.keys():
            residue_chunks[line[5]] = [line]
        else:
            residue_chunks[line[5]].append(line)

        if (
            line[3] in HALOGEN_ELECTRONEGATIVE.keys()
            and line[2] in HALOGEN_ELECTRONEGATIVE[line[3]]
        ):
            chain = line[4]
            if chain in chains:
                donor_atom_lines[chain].append(line)
        for key in HALOGENS.keys():
            # print(key)
            if str(key) in line[2]:
                halogen_atoms_lines.append(line)

    return halogen_atoms_lines, donor_atom_lines, residue_chunks


def make_pd_dataframe(final_interactions):
    df_dict = {
        "Focus Atom": [],
        "Focus Atom Res.Name": [],
        "Focus Atom Res.Number": [],
        "Focus Atom Chain": [],
        "Res. Name1": [],
        "Res. Number 1": [],
        "Chain 1": [],
        "Res. Name2": [],
        "Res. Number 2": [],
        "Chain 2": [],
        "Distance (Å)": [],
    }

    for instance in final_interactions.keys():
        df_dict["Focus Atom"].append(instance.split()[2])
        df_dict["Focus Atom Res.Name"].append(instance.split()[3])
        df_dict["Focus Atom Res.Number"].append(instance.split()[1])
        df_dict["Focus Atom Chain"].append(instance.split()[4])
        df_dict["Res. Name1"].append(final_interactions[instance][0][0][3])
        df_dict["Res. Number 1"].append(final_interactions[instance][0][0][1])
        df_dict["Chain 1"].append(final_interactions[instance][0][0][4])
        df_dict["Res. Name2"].append(final_interactions[instance][0][1][3])
        df_dict["Res. Number 2"].append(final_interactions[instance][0][1][1])
        df_dict["Chain 2"].append(final_interactions[instance][0][1][4])
        df_dict["Distance (Å)"].append(
            (final_interactions[instance][1][0], final_interactions[instance][1][1])
        )

    df = pd.DataFrame(df_dict)

    return df


def check_intercation_halogen(pdb_file, chains):
    halogen_atoms_lines, donor_atom_lines, residue_chunks = get_focussed_lines(
        pdb_file, chains
    )
    interactions = {}
    for halogen_line in halogen_atoms_lines:
        halogen_line1 = " ".join(str(e) for e in halogen_line)
        for chain in donor_atom_lines.keys():
            for donor_line in donor_atom_lines[chain]:
                dist = get_distance(halogen_line, donor_line)
                try:
                    if dist <= RADII[halogen_line[2]] + RADII[donor_line[2][0]]:
                        if halogen_line1 not in interactions.keys():
                            interactions[halogen_line1] = ([donor_line], [1, [dist]])
                        else:
                            interactions[halogen_line1][0].append(donor_line)
                            interactions[halogen_line1][1][0] += 1
                            interactions[halogen_line1][1][1].append(dist)
                except:
                    pass

    return interactions, residue_chunks


def get_X_Y_lines_from_residue(line, residue_chunk):
    # line_Y = [line]
    my_line = []
    dist = 100
    if type(line) is str:
        line = line.split(" ")
    for line1 in residue_chunk:
        if (
            line1 != line
            and get_distance(line, line1) < dist
            and get_distance(line, line1) != 0
        ):
            dist = get_distance(line, line1)
            my_line = line1
    return my_line


def angle_calculation(a, b, c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def calculate_angles(interactions: dict, residue_chunks: dict):

    for halogen_line in interactions.keys():
        halogen_Y_line = get_X_Y_lines_from_residue(
            halogen_line, residue_chunks[halogen_line.split()[5]]
        )
        for acceptor_line in interactions[halogen_line][0]:
            acceptor_Y_line = get_X_Y_lines_from_residue(
                acceptor_line, residue_chunks[acceptor_line[5]]
            )
            halogen_line1 = halogen_line.split(" ")

            point_halo = [
                float(halogen_line1[6]),
                float(halogen_line1[7]),
                float(halogen_line1[8]),
            ]
            point_acceptor_O = [
                float(acceptor_line[6]),
                float(acceptor_line[7]),
                float(acceptor_line[8]),
            ]
            point_halo_C = [
                float(halogen_Y_line[6]),
                float(halogen_Y_line[7]),
                float(halogen_Y_line[8]),
            ]

            theta1 = angle_calculation(point_halo_C, point_halo, point_acceptor_O)

            point_acceptor_Y = [
                float(acceptor_Y_line[6]),
                float(acceptor_Y_line[7]),
                float(acceptor_Y_line[8]),
            ]

            theta2 = angle_calculation(point_acceptor_Y, point_acceptor_O, point_halo)
            print(theta1, theta2)
            print(halogen_Y_line)

            if theta1 >= 120 and theta2 >= 110:
                interactions[halogen_line][1].append([theta1, theta2])

    return interactions


# interactions, residue_chunks = check_intercation_halogen("4i29.pdb", ["A", "B"])
# finalinteraction = calculate_angles(interactions, residue_chunks)

# for key in finalinteraction.keys():
#     print(key)
#     print(finalinteraction[key][0][0])
#     print(finalinteraction[key][1])
#     print()
