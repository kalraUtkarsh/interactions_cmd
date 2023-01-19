"""This file hosts the frontend of the Streamlit and all the workflow 
pipeline is set in this file"""
import os
import sys
import jinja2

import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
from matplotlib import pyplot as plt

import ressidues as RESIDUES
from inter_non_covalent_with_all_possibilities import make_files_for_all_rings
from segregation_all_possibilites import creating_sub_pdbs
from orient_non_covalent_refactor import read_sub_pdb
from filter_refactor import filtering
from Plotting import plotting
from water_mediated import check_intercation_water
from metal_mediated import check_intercation_metal
# from halogen_interaction import check_intercation_halogen
# from pi_pi_interactions import check_intercations_pi_pi
# from small_pdb_maker import make_small_pdb
# from cif2pdb import get_pdb_file


header = st.container()
results = st.container()

with header:
    st.image('coco.png')
st.sidebar.title("Specify your Parameters")

pdb_file = st.sidebar.file_uploader("Upload a pdb file")

chain1 = st.sidebar.text_area('Name of chain 1')
chain2 = st.sidebar.text_area('Name of chain 2')

threshold = st.sidebar.text_area('Threshold Value')

submit = send = st.sidebar.button("SUBMIT")

rings = RESIDUES.all_ring_names

types = ['cation_pi','anion_pi','sulfur_pi', 'lone_pair_pi'] #The types of interactions that are to be done with the same set of programmes

# print(pdb_file)
if pdb_file:
    fl = open(pdb_file.name, 'wb')
    fl.write(pdb_file.read())

    pdb_file = pdb_file.name


if submit:
    # if pdb_file.endswith('.cif'):
    #     get_pdb_file(pdb_file)
    #     pdb_file = pdb_file.replace('.cif','.pdb')
    placeholder = st.empty()
    placeholder.markdown("![Alt Text](https://media.giphy.com/media/pFZTlrO0MV6LoWSDXd/giphy.gif)")
    # make_small_pdb(pdb_file, [chain1, chain2], int(threshold))
    # pdb_file = f"{pdb_file}_small.pdb"
    output_file_name = f'{os.path.split(pdb_file)[1]}_feasible_ring_interactions.txt'
    
    chains = [chain1,chain2]

    new_dir = f'{os.path.split(pdb_file)[1]}_sub_pdb'
    os.makedirs(new_dir,exist_ok=True)
    os.makedirs(f'{new_dir}/Reoriented',exist_ok=True)
    os.makedirs(f'{new_dir}/Reoriented/xyz',exist_ok=True)
    os.makedirs(f'{new_dir}/Reoriented/data',exist_ok=True)
    os.makedirs(f'{new_dir}/Reoriented/filtered',exist_ok=True)
    print(output_file_name, pdb_file)

    # the specifications for the plotting 
    width = st.sidebar.slider("plot width", 1, 25, 3)
    height = st.sidebar.slider("plot height", 1, 25, 7)
    fig, ax = plt.subplots(figsize=(width, height))

    for type_ in types:
        output_file_with_interaction = open(output_file_name, 'a')
        print(pdb_file, output_file_with_interaction, type_)
        make_files_for_all_rings(rings,pdb_file, chains, output_file_with_interaction, type_)
    creating_sub_pdbs(output_file_name,pdb_file,new_dir)
    read_sub_pdb(new_dir)
    for type_ in types:   
        filtering(new_dir,type_)
    print('Done Cations and stuff')

    df_water = check_intercation_water(pdb_file, [chain1, chain2])
    print('Done Water')
    df_metal = check_intercation_metal(pdb_file, [chain1, chain2])
    print('Done Metal')
    # df_pi_pi = check_intercations_pi_pi(f'{new_dir}/Reoriented')
    # print('Done Pi-pi')

    extra_interactions = ['Water Mediated Interaction', 'Metal Mediated Interactions'] # Add other interactions to it when needed
    extra_df_s = [df_water, df_metal]

    print('Onto Plotting')
    
    """This section takes care of calling and processing of Plotting fucntions"""
    df, rest_dfs, plt1 = plotting(pdb_file, chains, threshold,f'{new_dir}/Reoriented/filtered', ax, fig, extra_df_s, extra_interactions)
    placeholder.empty()
    st.subheader(f'Plot of contacts with threshold {threshold}')
    st.pyplot(plt1)
    st.subheader('Table of contacts')
    st.dataframe(df)
    for rdf in list(rest_dfs.keys()):
        st.subheader(f"Table of interactions of type {rdf}")
        st.dataframe(rest_dfs[rdf])
    
    for extra, label_ in zip(extra_df_s, extra_interactions):
        print(type(extra))
        if extra is not None:
            st.subheader(f"Table of {label_}")
            st.dataframe(extra)
