import os

from read_pdb_regex import read_pdb

from Bio.PDB import PDBParser, PDBIO

io = PDBIO()


def get_all_chains(pdb_file):
    pdb = PDBParser().get_structure(pdb_file, pdb_file)
    chains = []
    for chain in pdb.get_chains():
        io.set_structure(chain)
        chains.append(chain.get_id())
        # io.save(pdb.get_id() + "_" + chain.get_id() + ".pdb")
    return chains
