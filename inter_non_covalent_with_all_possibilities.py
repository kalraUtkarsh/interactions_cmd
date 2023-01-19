import os
import sys

from read_pdb_regex import read_pdb
from read_pdb_regex import distance
import ressidues as RESIDUES


"""Function that searches for all possible interactions 
of a FOCUS_ATOM with a ring.

param rings: name of all the aromatic rings possible. TYPE: list 
param pdb_file: path of the pdb_file                  TYPE: string
"""


def make_files_for_all_rings(
    rings: list, pdb_file: str, chains: list, output_file_with_interaction, type_
):

    # print('WE')
    parsed, lines = read_pdb(f"./{pdb_file}")
    for ring in rings:

        ring_atoms = RESIDUES.residues[ring]

        aminos = list(getattr(RESIDUES, type_).keys())
        for amino in aminos:
            for focus_atom in getattr(RESIDUES, type_)[amino]:
                # print(focus_atom)
                a_data = []
                # print(focus_atom)
                for parsed_line in parsed:
                    # print(parsed_line[2], focus_atom, parsed_line[3], amino, parsed_line[4], chains)
                    # print(parsed_line)
                    if type_ == "halogen_pi":
                        if (
                            parsed_line[2] == focus_atom
                            and parsed_line[3] != amino
                            and parsed_line[4] in chains
                        ):
                            # print(parsed_line)
                            a_data.append(parsed_line)
                    else:
                        if (
                            parsed_line[2] == focus_atom
                            and parsed_line[3] == amino
                            and parsed_line[4] in chains
                        ):
                            # print(parsed_line)
                            a_data.append(parsed_line)
                # print(a_data)
                matched = matching_function(a_data, parsed, ring, ring_atoms, chains)
                # print(matched)
                # print('HERE')
                dictance_checking_function(
                    matched, os.path.split(pdb_file)[1], output_file_with_interaction
                )
    output_file_with_interaction.close()


"""Function that needs to be called from inside of the 
make_files_for_all_rings function. It checks for the pdb line that
has the matching ring and atoms we are looking for

param a_data: list of all the lines of the pdb file that has the FOCUS_ATOM
              with the same amino acid.
param parsed: all lines of the pdb file.
param ring: name of the current ring we are looking for.
param ring_atoms: list of all the atoms that aromatic ring posseses."""


def matching_function(a_data, parsed, ring, ring_atoms, chains):

    matched = {}

    for a in a_data:
        for i in parsed:

            if (
                i[3] == ring and i[2] in ring_atoms and i[4] in chains
            ):  # and i[4] != a[4]:
                # print(i)
                # key = a[11] + '_' + '_'.join(a[4:6]) + ' | ' + '_'.join(i[3:6])
                key = "_".join(a[3:6]) + " | " + "_".join(i[3:6])
                # print(key)
                if key not in matched.keys():
                    matched[key] = []
                    matched[key].append(a)
                if i not in matched[key]:
                    matched[key].append(i)

    return matched


"""checks if the dictance between the given atoms is less than 9.0
param matched: list of matched lines
param file: name of the pdb file
param fw: name of the file in which to be written"""


def dictance_checking_function(matched, file, fw):

    for key in matched.keys():
        to_check = matched[key]
        if len(to_check) == 7:
            # print(to_check)
            if (
                distance(to_check[0], to_check[1]) < 9.0
                and distance(to_check[0], to_check[2]) < 9.0
                and distance(to_check[0], to_check[3]) < 9.0
                and distance(to_check[0], to_check[4]) < 9.0
                and distance(to_check[0], to_check[5]) < 9.0
                and distance(to_check[0], to_check[6]) < 9.0
            ):
                # print(to_check)
                fw.write(file + " | " + key + "\n")
                # print(file + ' ' + key + '\n')
        elif len(to_check) == 10:
            if (
                distance(to_check[0], to_check[1]) < 9.0
                and distance(to_check[0], to_check[2]) < 9.0
                and distance(to_check[0], to_check[3]) < 9.0
                and distance(to_check[0], to_check[4]) < 9.0
                and distance(to_check[0], to_check[5]) < 9.0
                and distance(to_check[0], to_check[6]) < 9.0
                and distance(to_check[0], to_check[7]) < 9.0
                and distance(to_check[0], to_check[8]) < 9.0
                and distance(to_check[0], to_check[9]) < 9.0
            ):
                # print(to_check)
                fw.write(file + " | " + key + "\n")
                # print(file + ' ' + key + '\n')


# rings = RESIDUES.all_ring_names
# pdb_file = "1p54.pdb"
# chains = ["A", "B"]
# output_file_with_interaction = open("try1.txt", "w")
# type_ = "halogen_pi"
# make_files_for_all_rings(rings, pdb_file, chains, output_file_with_interaction, type_)
