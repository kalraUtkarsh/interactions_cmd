"""This file takes care of the Plotting. The Plotting is still not fully done
There can be some minor problems"""
# from cProfile import label
import sys
import os

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from read_pdb_regex import read_pdb

plt.rcParams["figure.figsize"] = (7,20)

def distance_table_creator(chains):
    distance_table = []
    distances = []
    threshold = 9
    for residue1 in chains[list(chains.keys())[0]]:
        for residue2 in chains[list(chains.keys())[1]]:
            x1 = float(residue1[6])
            y1 = float(residue1[7])
            z1 = float(residue1[8])
            x2 = float(residue2[6])
            y2 = float(residue2[7])
            z2 = float(residue2[8])

            a = np.array((x1, y1, z1))
            b = np.array((x2, y2, z2))
            dist = np.linalg.norm(a-b)
            distances.append(dist)
            if dist <= threshold:
                line = []
                # line.append(residue1[2])
                line.append(residue1[3])
                line.append(residue1[5])
                line.append(residue1[4])
                # line.append(residue2[2])
                line.append(residue2[3])
                line.append(residue2[5])
                line.append(residue2[4])
                line.append(dist)
                distance_table.append(line)
    # print(distance_table)
    return distance_table


def initial_filter_data_frame_creator(pdb_file, chains_, threshold):
    if chains_[0] != chains_[1]:
        chains = {f'{chains_[0]}':[],f'{chains_[1]}':[]}
        if pdb_file.endswith('.pdb'):
            fl = open(pdb_file,'r').readlines()

            for line in fl:
                line = line.split()
                # print(line)
                if 'ATOM' == line[0]:
                    # print(line)
                    if line[4] in chains_:
                        # print(line)
                        chain = line[4]
                        chains[chain].append(line)
        if pdb_file.endswith('.cif'):
            lines, emp = read_pdb(pdb_file)
            for line in lines:
                line.remove('')
                line.remove(' ')
                # print(line)
                if 'ATOM' == line[0]:
                    if line[4] in list(chains.keys()):
                        
                        chain = line[4]
                        chains[chain].append(line)
    else:
        chains = {f'{chains_[0]}':[],f'{chains_[1]}*':[]}
        if pdb_file.endswith('.pdb'):
            fl = open(pdb_file,'r').readlines()

            for line in fl:
                line = line.split()
                # print(line)
                if 'ATOM' == line[0]:
                    # print(line)
                    if line[4] in chains_:
                        # print(line)
                        chain = line[4]
                        chains[chain].append(line)
                        chains[f'{chain}*'].append(line)
        if pdb_file.endswith('.cif'):
            lines, emp = read_pdb(pdb_file)
            for line in lines:
                line.remove('')
                line.remove(' ')
                # print(line)
                if 'ATOM' == line[0]:
                    if line[4] in list(chains.keys()):
                        
                        chain = line[4]
                        chains[chain].append(line)
                        chains[f'{chain}*'].append(line)
    # print(chains)
    # print(chains[chains_[1]][0])
    print(chains[list(chains.keys())[0]][0][5])
    x_begin, x_end = int(chains[list(chains.keys())[0]][0][5]), int(chains[list(chains.keys())[0]][-1][5])
    y_begin, y_end = int(chains[list(chains.keys())[1]][0][5]), int(chains[list(chains.keys())[1]][-1][5])
    
    distance_table = distance_table_creator(chains)

    columns = ['Res. Name1', 'Res. Number 1', 'Chain 1','Res. Name2', 'Res. Number 2', 'Chain 2','Distance (Å)']
    df = pd.DataFrame(distance_table)
    df.columns = columns
    df = df.drop_duplicates(subset=['Res. Name1', 'Res. Number 1', 'Chain 1','Res. Name2', 'Res. Number 2', 'Chain 2'], keep="last")

    return df, x_begin, x_end, y_begin, y_end
    


def plotting(pdb_file, chains, threshold, filtered_dir, ax, fig, extra_df_s, extra_interactions):

    df, x_begin, x_end, y_begin, y_end = initial_filter_data_frame_creator(pdb_file, chains, threshold)

    df.to_csv('bjla.csv')
    x = df['Res. Number 1'].tolist()
    y = df['Res. Number 2'].tolist()
    # print(x)
    # print(y)
    X = x
    Y = y
    x = []
    y = [] 
    for i,j in zip(X,Y):
        try:
            x.append(int(i))
            y.append(int(j))
        except:
            pass


    rest_plots = {}
    rest_dfs = {}
    for index,instance in enumerate(os.listdir(filtered_dir)):
        rest_plots[instance] = {'x':[],'y':[]}
        rest_dfs[instance] = {'Res. Name1':[], 'Res. Number 1':[], 'Chain 1':[],'Res. Name2':[], 'Res. Number 2':[], 'Chain 2':[],'Distance (Å)':[]}
        for index1, instance1 in enumerate(os.listdir(f'{filtered_dir}/{instance}')):
            chain1 = instance1.split('-')[0].split('_')[1]
            x1  = int(instance1.split('-')[0].split('_')[2])
            y1 = int(instance1.split('-')[1].split('_')[2])
            n1 = instance1.split('-')[0].split('_')[0]
            n2 = instance1.split('-')[1].split('_')[0]
            c1 = instance1.split('-')[0].split('_')[1]
            c2 = instance1.split('-')[1].split('_')[1]
            d = np.nan
            if chain1 == chains[0]:
                rest_plots[instance]['x'].append(x1)
                rest_plots[instance]['y'].append(y1)
                rest_dfs[instance]['Res. Name1'].append(n1)
                rest_dfs[instance]['Res. Name2'].append(n2)
                rest_dfs[instance]['Res. Number 1'].append(x1)
                rest_dfs[instance]['Res. Number 2'].append(y1)
                rest_dfs[instance]['Chain 1'].append(c1)
                rest_dfs[instance]['Chain 2'].append(c2)
                rest_dfs[instance]['Distance (Å)'].append(d)
            else:
                rest_plots[instance]['x'].append(y1)
                rest_plots[instance]['y'].append(x1)
                rest_dfs[instance]['Res. Name1'].append(n2)
                rest_dfs[instance]['Res. Name2'].append(n1)
                rest_dfs[instance]['Res. Number 1'].append(y1)
                rest_dfs[instance]['Res. Number 2'].append(x1)
                rest_dfs[instance]['Chain 1'].append(c2)
                rest_dfs[instance]['Chain 2'].append(c1)
                rest_dfs[instance]['Distance (Å)'].append(d)
        
        rest_dfs[instance] = pd.DataFrame(rest_dfs[instance])


    ax.scatter(x, y,  color='black')
    for graph in list(rest_plots.keys()):
        ax.plot(rest_plots[graph]['x'],rest_plots[graph]['y'], '^',label = graph)
    
    for extra, label_ in zip(extra_df_s, extra_interactions):
        print(type(extra))
        if extra is not None:
            print(label_)
            ax.plot([int(x) for x in extra['Res. Number 1'].tolist()],[int(x) for x in extra['Res. Number 2'].tolist()],'o',label = label_)

    # ax.xlim([x_begin,x_end])
    # ax.ylim([y_begin,y_end])
    ax.axis(xmin=x_begin, ymin=y_begin)
    ax.axis(xmax=x_end, ymax=y_end)
    # ax.ylabel(chains[1])
    # plt.xlabel(chains[0])
    fig.supxlabel(chains[0])
    fig.supylabel(chains[1])
    ax.legend()
    # plt.show()
    print(rest_plots)
    return df, rest_dfs, fig


def try_():
    x = [i for i in range(10)]
    plt.plot(x,x)
    plt.show()


# initial_filter_data_frame_creator('1pue.pdb', ['E','E'], 6)
