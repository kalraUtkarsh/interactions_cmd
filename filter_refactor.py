import os 

from shutil import copyfile, copy

import ressidues as RESIDUES


def filtering(sub_pdb_dir, type_):

    for data_file in os.listdir(f'{sub_pdb_dir}/Reoriented/data'):
        ring = data_file.split('_')[0]
        amino = data_file.split('-')[1].split('_')[0]
        a = RESIDUES.THRESHOLDS[RESIDUES.axis_definition_dictionary[ring][2]][0]
        b = RESIDUES.THRESHOLDS[RESIDUES.axis_definition_dictionary[ring][2]][1]
        if amino in list(getattr(RESIDUES,type_).keys()):
            molecules = getattr(RESIDUES,type_)[amino]
            # print(molecules, amino)
            data_file = f'{sub_pdb_dir}/Reoriented/data/{data_file}' 
            for molecule in molecules:
                with open(data_file) as fh:
                    for line in fh.readlines():
                        if line.startswith(molecule):
                            # print(line)
                            atom = line.split(' | ')[0]
                            x = float(line.split(' | ')[1])
                            y = float(line.split(' | ')[2])
                            z = float(line.split(' | ')[3])
                            dist = x*x/a + y*y/b
                            # print(dist,abs(z))
                            if dist <= 1 and abs(z) <= 6: # here is the rise parameter
                                # print('yes')
                                # print(data_file, '\n', atom, x, y, z, x*x + y*y)
                                src_dir = f'{sub_pdb_dir}/Reoriented/xyz/'
                                src = src_dir + data_file.split('/')[-1].split('-re.data')[0] + '-reoriented.xyz'
                                if not os.path.exists(os.getcwd() +'/'+ sub_pdb_dir +f'/Reoriented/filtered/{type_}'):
                                    os.mkdir(os.getcwd() +'/'+ sub_pdb_dir +f'/Reoriented/filtered/{type_}')
                                
                                dst = os.getcwd() +'/'+ sub_pdb_dir +f'/Reoriented/filtered/{type_}'
                                copy(src, dst)

# filtering("4k50.pdb_small.pdb_sub_pdb","lone_pair_pi")