"""This file hosts the frontend of the Streamlit and all the workflow 
pipeline is set in this file"""
import os
import sys
import jinja2

import pandas as pd
from matplotlib import pyplot as plt
import fire

import ressidues as RESIDUES
from inter_non_covalent_with_all_possibilities import make_files_for_all_rings
from segregation_all_possibilites import creating_sub_pdbs
from orient_non_covalent_refactor import read_sub_pdb
from filter_refactor import filtering
from Plotting import plotting
from water_mediated import check_intercation_water
from metal_mediated import check_intercation_metal
from get_all_chains import get_all_chains

# from halogen_interaction import check_intercation_halogen
# from pi_pi_interactions import check_intercations_pi_pi
# from small_pdb_maker import make_small_pdb
# from cif2pdb import get_pdb_file


rings = RESIDUES.all_ring_names

types = [
    "cation_pi",
    "anion_pi",
    "sulfur_pi",
    "lone_pair_pi",
]  # The types of interactions that are to be done with the same set of programmes

# print(pdb_file)
def run_programme(pdb_file: str, mode: str, chain1=None, chain2=None, threshold=6):
    all_chains = []
    if mode == "all":
        all_chains = get_all_chains(pdb_file)
    # if pdb_file.endswith('.cif'):
    #     get_pdb_file(pdb_file)
    #     pdb_file = pdb_file.replace('.cif','.pdb')
    # make_small_pdb(pdb_file, [chain1, chain2], int(threshold))
    # pdb_file = f"{pdb_file}_small.pdb"
    output_file_name = f"{os.path.split(pdb_file)[1]}_feasible_ring_interactions.txt"

    chains = [chain1, chain2]

    new_dir = f"{os.path.split(pdb_file)[1]}_sub_pdb"
    os.makedirs(new_dir, exist_ok=True)
    os.makedirs(f"{new_dir}/Reoriented", exist_ok=True)
    os.makedirs(f"{new_dir}/Reoriented/xyz", exist_ok=True)
    os.makedirs(f"{new_dir}/Reoriented/data", exist_ok=True)
    os.makedirs(f"{new_dir}/Reoriented/filtered", exist_ok=True)
    print(output_file_name, pdb_file)

    # the specifications for the plotting
    fig, ax = plt.subplots(figsize=(5, 10))

    for type_ in types:
        output_file_with_interaction = open(output_file_name, "a")
        print(pdb_file, output_file_with_interaction, type_)
        make_files_for_all_rings(
            rings, pdb_file, chains, output_file_with_interaction, type_
        )
    creating_sub_pdbs(output_file_name, pdb_file, new_dir)
    read_sub_pdb(new_dir)
    for type_ in types:
        filtering(new_dir, type_)
    print("Done Cations and stuff")

    df_water = check_intercation_water(pdb_file, [chain1, chain2])
    print("Done Water")
    df_metal = check_intercation_metal(pdb_file, [chain1, chain2])
    print("Done Metal")
    # df_pi_pi = check_intercations_pi_pi(f'{new_dir}/Reoriented')
    # print('Done Pi-pi')

    extra_interactions = [
        "Water Mediated Interaction",
        "Metal Mediated Interactions",
    ]  # Add other interactions to it when needed
    extra_df_s = [df_water, df_metal]

    print("Onto Plotting")

    """This section takes care of calling and processing of Plotting fucntions"""
    df, rest_dfs, plt1 = plotting(
        pdb_file,
        chains,
        threshold,
        f"{new_dir}/Reoriented/filtered",
        ax,
        fig,
        extra_df_s,
        extra_interactions,
    )
    plt1.savefig("interactions.png")
    df.to_excel("table_of_contents.xlsx")
    for rdf in list(rest_dfs.keys()):
        rest_dfs[rdf].to_excel(f"interactions_of_type_{rdf}.xlsx")

    for extra, label_ in zip(extra_df_s, extra_interactions):
        print(type(extra))
        if extra is not None:
            extra.to_excel(f"Table_of_{label_}.xlsx")


if __name__ == "__main__":
    fire.Fire(run_programme)
