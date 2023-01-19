# from cProfile import label
# import sys
# import os

import pandas as pd
import numpy as np
# from matplotlib import pyplot as plt
from read_pdb_regex import read_pdb


def make_small_pdb(pdb_file, chains_, threshold):
    chains = {chains_[0]:[],chains_[1]:[]}
    if pdb_file.endswith('.pdb'):
        fl = open(f"{pdb_file}_small.pdb",'w')
        parsed, lines = read_pdb(pdb_file)
        uninterupted = []
        for line1, line2 in zip(parsed, lines):
            written_check = 0
            if line1[4] in chains_:
                fl.write(f'{line2}\n')
                # for line3, line4 in zip(parsed, lines):
                #     if line3[4] in chains_:
                #         x1 = float(line1[6])
                #         y1 = float(line1[7])
                #         z1 = float(line1[8])
                #         x2 = float(line3[6])
                #         y2 = float(line3[7])
                #         z2 = float(line3[8])
                #         a = np.array((x1, y1, z1))
                #         b = np.array((x2, y2, z2))
                #         dist = np.linalg.norm(a-b)

                #         if dist <= threshold:
                #             if line1[-1] != 0:
                #                 line2 = f"{line2}\n"
                #                 fl.write(line2)
                #                 line1.append(0)
                #             if line3[-1] != 0:
                #                 line4 = f"{line4}\n"
                #                 fl.write(line4)
                #                 line3.append(0)
        fl.close()
                            
make_small_pdb('3ts0.pdb', ['U','A'], 6)