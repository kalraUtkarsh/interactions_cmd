# Theresholds to be 4 ang, 45 degrees for theta. x_ring1_max < x_ring2_min
import os
import math

import numpy as np

import ressidues as RESIDUES

rings = RESIDUES.residues.keys()
rings_with_atoms = RESIDUES.residues

"""General functions of this file includes taking in all the 
sub pdbs generated during the cation pi and those kind of interaction 
process and then checking the Pi-Pi conditions on that small files."""

def matching_rings(dir_):
    """This function takes in the directory and creates a dictionary 
    of the rings, structure of this dictionary is given in the 
    description of next function"""

    matched = {}
    for file1, file2 in zip(os.listdir(f"{dir_}"),os.listdir(f"{dir_}")):
        if file1.endswith('.pdb'):
            # print(file1)
            ring1 = file1.split('_')[0]
            ring2 = file1.split('-')[1].split('_')[0]

            if ring1 in rings and ring2 in rings:
                ring_coordinates = {ring1:{}, ring2:{}}
                
                file1_open = open(f"{dir_}/{file1}",'r').readlines()
                file2_open = open(f"{dir_}/{file2}",'r').readlines()

                for index, line in enumerate(file2_open):
                    atom = line.split()[2].strip()
                    Ring = line.split()[3]
                    coords = file1_open[index].split()
                    if atom in rings_with_atoms[ring1] and Ring == ring1:
                        # print(atom)
                        ring_coordinates[ring1][atom]=[float(coords[6].strip()),float(coords[7].strip()),float(coords[8].strip())]
                        # print(ring_coordinates[ring1][atom])

                    if atom in rings_with_atoms[ring2] and Ring == ring2:
                        ring_coordinates[ring2][atom]=[float(coords[6].strip()),float(coords[7].strip()),float(coords[8].strip())]

                matched[file1] = ring_coordinates
    
    return matched

# def create_dataframe(final_files):



def checking_interactions_pi_pi(matched, threshold):
    """This funtion takes in the dictionary created in the previous 
    step and then calculates the angles between the two rings.
    Namely theta and gama"""

    """Param matched: dictionary of rings
       Param threshold: threshold of distance between the centroids of 2 rings
       Return fie_names: list of lists [['name of the file', Theta, gama, centroind_distance]]"""

    """matched={name of the file:
                                {ring1:
                                        {atom1:[x,y,z],
                                        atom2:[x,y,z]}
                                },
                                ring2:
                                        {atom1:[x,y,z],
                                        atom2:[x,y,z]}
                }"""

    centroids = {}
    final_files = []
    for file_ in matched.keys():
        # print(file_)
        # print(matched[file_].keys())
        try:
            centroids[file_] = {}
            for ring in matched[file_].keys():
                # print(ring)
                
                centroids[file_][ring] = {'x':0,'y':0,'z':0}
                for atom in matched[file_][ring].keys():
                    centroids[file_][ring]['x'] += matched[file_][ring][atom][0]
                    centroids[file_][ring]['y'] += matched[file_][ring][atom][1]
                    centroids[file_][ring]['z'] += matched[file_][ring][atom][2]
                centroids[file_][ring]['x'] /= len(matched[file_][ring].keys())
                centroids[file_][ring]['y'] /= len(matched[file_][ring].keys())
                centroids[file_][ring]['z'] /= len(matched[file_][ring].keys())

            rings = list(centroids[file_].keys())
            a = np.array([centroids[file_][rings[0]]['x'],centroids[file_][rings[0]]['y'],centroids[file_][rings[0]]['z']])
            b = np.array([centroids[file_][rings[1]]['x'],centroids[file_][rings[1]]['y'],centroids[file_][rings[1]]['z']])
            r_cen = np.linalg.norm(a-b)
            print(r_cen)
            if r_cen <= threshold:
                # print()
                # print(file_)
                # print(f"Distance between {rings[0]} and {rings[1]} is {r_cen}")
                atoms_ring1 = list(matched[file_][rings[0]].keys())
                p1 = a
                p2 = np.array(matched[file_][rings[0]][atoms_ring1[0]])
                p3 = np.array(matched[file_][rings[0]][atoms_ring1[1]])

                v1 = p3 - p1
                v2 = p2 - p1

                # the cross product is a vector normal to the plane
                normal_vector = np.cross(v1, v2)
                a1, b1, c1 = normal_vector
                centroids_vector = a - b
                dot_prod = np.linalg.norm(np.dot(normal_vector,centroids_vector))
                cos_theta = dot_prod/(r_cen*math.sqrt(a1*a1+b1*b1+c1*c1)) 
                theta = np.degrees(np.arccos(cos_theta))
                # print(theta, 'THETA')


                atoms_ring2 = list(matched[file_][rings[1]].keys())
                p1 = b
                p2 = np.array(matched[file_][rings[1]][atoms_ring2[0]])
                p3 = np.array(matched[file_][rings[1]][atoms_ring2[1]])

                v1 = p3 - p1
                v2 = p2 - p1

                # the cross product is a vector normal to the plane
                normal_vector_ = np.cross(v1, v2)
                a1, b1, c1 = normal_vector_
                centroids_vector_ = a - b
                dot_prod = np.linalg.norm(np.dot(normal_vector_,normal_vector))
                cos_theta_ = dot_prod/(r_cen*math.sqrt(a1*a1+b1*b1+c1*c1)) 
                gama = np.degrees(np.arccos(cos_theta_))
                final_files.append([file_,r_cen,theta, gama])
        except:
            pass
        
    return final_files





# print(matching_rings("/mnt/d/ProtTest2/cocomaps/4k50.pdb_smalll.pdb_sub_pdb/Reoriented"))
sd = open('pipi.txt','a')                   
for res in checking_interactions_pi_pi(matching_rings("D:\ProtTest2\cocomaps\\1pue.pdb_sub_pdb"),6):
    print(res)
    sd.write(f'{res}\n')
    





