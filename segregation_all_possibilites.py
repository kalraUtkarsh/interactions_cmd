import os

from collections import Counter

import ressidues as RESIDUES
from read_pdb_regex import read_pdb
from read_pdb_regex import distance

# pdb_file_parsed = read_pdb(pdb_file)
# output_file_with_interaction = '1awc.pdb_feasible_ring_interactions.txt'

"""Reads the file created in the previous step. segregates 
and returns the information in a structured format.

param lines: all line of the output file of last step
return all_occourances: format: list of lists of unique occourances
return residue_numbers: format: Counter dict of residue numbers"""


def read_output_file_lines(lines):
    all_occourances = []
    all_occourances_crude = []
    residues_numbers = []
    for line in lines:
        line = line.strip("\n")
        if line not in all_occourances_crude:
            all_occourances_crude.append(line)
            line = line.split("|")
            amino = line[1].split("_")[0].strip()
            chain1 = line[1].split("_")[1].strip()
            residue_number_1 = int(line[1].split("_")[2].strip())
            ring_name = line[2].split("_")[0].strip()
            chain2 = line[2].split("_")[1].strip()
            residue_number_2 = int(line[2].split("_")[2].strip())
            all_occourances.append(
                [amino, chain1, residue_number_1, ring_name, chain2, residue_number_2]
            )
            residues_numbers.append(residue_number_1)
            residues_numbers.append(residue_number_2)
    residues_numbers.sort()
    residues_numbers = dict(Counter(residues_numbers))
    return all_occourances, residues_numbers


"""creates sub pdbs of the indentified possible interactoins.

param output_file_with_interaction: file name of the output file"""


def creating_sub_pdbs(output_file_with_interaction, pdb_file, sub_pdb_dir):
    possibilities = open(output_file_with_interaction, "r").readlines()
    list_of_residues_to_read_from_pdb, residues_numbers = read_output_file_lines(
        possibilities
    )
    parsed_lines, bla = read_pdb(pdb_file)
    for index1, residue in enumerate(list_of_residues_to_read_from_pdb):
        file_name = f"{sub_pdb_dir}/{residue[3]}_{residue[4]}_{residue[5]}-{residue[0]}_{residue[1]}_{residue[2]}-{os.path.basename(pdb_file)}"
        fl = open(file_name, "w")
        res1, amino1, chain1, res2, amino2, chain2 = "", "", "", "", "", ""

        if min(residue[2], residue[5]) == residue[5]:
            res1 = residue[5]
            amino1 = residue[3]
            chain1 = residue[4]
            res2 = residue[2]
            amino2 = residue[0]
            chain2 = residue[1]
        else:
            res1 = residue[2]
            amino1 = residue[0]
            chain1 = residue[1]
            res2 = residue[5]
            amino2 = residue[3]
            chain2 = residue[4]

        for index2, line in enumerate(parsed_lines):
            if line[5] == str(res1) and line[3] == amino1 and line[4] == chain1:
                printing_line = ""
                for i in line:
                    printing_line += str(i) + "\t"
                # printing_line = printing_line.strip()+'\n'
                printing_line = bla[index2] + "\n"
                fl.write(printing_line)
            if line[5] == str(res2) and line[3] == amino2 and line[4] == chain2:
                printing_line = ""
                for i in line:
                    printing_line += str(i) + "\t"
                # printing_line = printing_line.strip()+'\n'
                printing_line = bla[index2] + "\n"
                fl.write(printing_line)


# creating_sub_pdbs("try1.txt", "./1p54.pdb", "./1p54_halo")
