import os
import numpy as np
import math

import ressidues as RESIDUES
from read_pdb_regex import read_pdb


def magnitude(arr, d1, d2):
    result = math.sqrt(arr[d1] * arr[d1] + arr[d2] * arr[d2])
    return result


def read_sub_pdb(dir):
    atom = 2

    for sub_pdb in os.listdir(dir):
        if sub_pdb.endswith(".pdb"):
            centroid = np.array([0.0, 0.0, 0.0])
            xaxis = np.array([0.0, 0.0, 0.0])
            yaxis = np.array([0.0, 0.0, 0.0])
            coord = []

            ring = sub_pdb.split("_")[0]
            # amino = sub_pdb.split('-')[1].split('_')[0]

            x_axis_atoms = RESIDUES.axis_definition_dictionary[ring][0]
            y_axis_atoms = RESIDUES.axis_definition_dictionary[ring][1]
            div = RESIDUES.axis_definition_dictionary[ring][2]
            ring_atoms = RESIDUES.residues[ring]

            sub_pdb_data, full_line = read_pdb(f"{dir}/{sub_pdb}")

            for row in sub_pdb_data:
                tmp = np.array(row[6:9])
                coord.append(tmp)
                if row[2] in ring_atoms:
                    # print(row)
                    centroid += tmp
                    if row[atom] in x_axis_atoms:
                        xaxis += tmp
                    if row[atom] in y_axis_atoms:
                        yaxis += tmp

            centroid /= div
            xaxis /= len(x_axis_atoms)
            coord = np.array(coord)

            for store in [coord, xaxis, yaxis, centroid]:
                store -= centroid

            zaxis = np.cross(xaxis, yaxis)
            yaxis = np.cross(zaxis, xaxis)

            x = 0
            y = 1
            z = 2

            # First rotation along Z axis, compute the angle and transform
            z_mag = magnitude(xaxis, x, y)
            z_cos = xaxis[x] / z_mag
            z_sin = xaxis[y] / z_mag
            # print(z_mag, z_cos, z_sin)

            for i in coord:
                x_tmp = i[x] * z_cos + i[y] * z_sin
                y_tmp = -i[x] * z_sin + i[y] * z_cos
                i[x] = x_tmp
                i[y] = y_tmp

            for i in [xaxis, yaxis]:
                x_tmp = i[x] * z_cos + i[y] * z_sin
                y_tmp = -i[x] * z_sin + i[y] * z_cos
                i[x] = x_tmp
                i[y] = y_tmp

            # print('After Z')
            # print(xaxis)
            # print(yaxis)

            # Second rotation along Y axis
            y_mag = magnitude(xaxis, x, z)
            y_cos = xaxis[x] / y_mag
            y_sin = xaxis[z] / y_mag
            # print(y_mag, y_cos, y_sin)
            for i in coord:
                x_tmp = i[x] * y_cos + i[z] * y_sin
                z_tmp = -i[x] * y_sin + i[z] * y_cos
                i[x] = x_tmp
                i[z] = z_tmp

            for i in [xaxis, yaxis]:
                x_tmp = i[x] * y_cos + i[z] * y_sin
                z_tmp = -i[x] * y_sin + i[z] * y_cos
                i[x] = x_tmp
                i[z] = z_tmp

            # print('After Y')
            # print(xaxis)
            # print(yaxis)

            # Third rotation along X axis
            x_mag = magnitude(yaxis, y, z)
            x_cos = yaxis[y] / x_mag
            x_sin = yaxis[z] / x_mag
            # print(x_mag, x_cos, x_sin)
            for i in coord:
                y_tmp = i[y] * x_cos + i[z] * x_sin
                z_tmp = -i[y] * x_sin + i[z] * x_cos
                i[y] = y_tmp
                i[z] = z_tmp

            for i in [xaxis, yaxis]:
                y_tmp = i[y] * x_cos + i[z] * x_sin
                z_tmp = -i[y] * x_sin + i[z] * x_cos
                i[y] = y_tmp
                i[z] = z_tmp

            create_reoriented_sub(dir, sub_pdb, coord)


def create_reoriented_sub(dir, file, coordinates):
    refile = file.split(".")[0] + "-reoriented.xyz"
    refile = f"{dir}/Reoriented/xyz/{refile}"
    dt = file.split(".")[0] + "-re.data"
    dt = f"{dir}/Reoriented/data/{dt}"
    data, somthing = read_pdb(f"{dir}/{file}")
    with open(refile, "w") as fout:
        with open(dt, "w") as fdt:
            counter = 0
            fout.write(str(len(data)) + "\n")
            fdt.write(str(len(data)) + "\n")
            fout.write(file + "\n")
            fdt.write(file + "\n")
            for row in data:
                for i in range(len(coordinates[counter])):
                    coordinates[counter][i] = (
                        math.ceil(coordinates[counter][i] * 100000) / 100000.0
                    )
                row[7:10] = coordinates[counter]
                row = [str(i) for i in row]
                to_file = (
                    row[2][0]
                    + "     "
                    + row[7].rjust(8)
                    + "  "
                    + row[8].rjust(8)
                    + "  "
                    + row[9].rjust(8)
                )
                fout.write(to_file + "\n")
                to_dt = row[2] + " | " + row[7] + " | " + row[8] + " | " + row[9]
                fdt.write(to_dt + "\n")
                counter += 1


read_sub_pdb("1p54_halo")
