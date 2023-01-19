import os 

from shutil import copyfile, copy

import ressidues as RESIDUES


def filtering(sub_pdb_dir, type_):

    for data_file in os.listdir(f'{sub_pdb_dir}/Reoriented/data'):
        ring = data_file.split('_')[0]    # A
        molecule = data_file.split('-')[1].split('_')[0] # HIS
        a = RESIDUES.THRESHOLDS[RESIDUES.axis_definition_dictionary[ring][2]][0]
        b = RESIDUES.THRESHOLDS[RESIDUES.axis_definition_dictionary[ring][2]][1]
        X_atoms = RESIDUES.axis_definition_dictionary[ring][0] # ['N1','C2']
        Y_atom = RESIDUES.axis_definition_dictionary[ring][1] # C6
        molecule_atoms = RESIDUES.getattr(RESIDUES,type_)[molecule] # ['O']

        X_atoms_coordinates = []
        Y_atoms_coordinates= []
        fh = open(f'{sub_pdb_dir}/Reoriented/data/{data_file}', 'r')
        for line in fh.readlines():
            potential_atom = line.split()[0]
            if potential_atom in X_atoms:
                x = float(line.split(' | ')[1])
                y = float(line.split(' | ')[2])
                z = float(line.split(' | ')[3])

                X_atoms_coordinates.append([x,y,z])

            if potential_atom == Y_atom:
                x = float(line.split(' | ')[1])
                y = float(line.split(' | ')[2])
                z = float(line.split(' | ')[3])

                Y_atoms_coordinates.append([x,y,z])
        
        Final_coordinates_X = [0,0,0]
        for x in X_atoms_coordinates:
            Final_coordinates_X[0] += Final_coordinates_X[0] + x[0]
            Final_coordinates_X[1] += Final_coordinates_X[1] + x[1]
            Final_coordinates_X[2] += Final_coordinates_X[2] + x[2]
        Final_coordinates_X = [x/len(X_atoms_coordinates) for x in Final_coordinates_X]

        Final_coordinates_Y = [0,0,0]
        for x in X_atoms_coordinates:
            Final_coordinates_Y[0] += Final_coordinates_Y[0] + x[0]
            Final_coordinates_Y[1] += Final_coordinates_Y[1] + x[1]
            Final_coordinates_Y[2] += Final_coordinates_Y[2] + x[2]
        Final_coordinates_Y = [x/len(X_atoms_coordinates) for x in Final_coordinates_Y]

        



        