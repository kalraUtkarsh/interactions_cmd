import ressidues as RESIDUES
from read_pdb_regex import read_pdb


def get_residue_chunk(pdb_file, res_name):
    parsed, lines = read_pdb(pdb_file)
    residue_chunks = {}

    for line in parsed:
        if line[3] not in residue_chunks.keys() and line[3] == res_name:
            residue_chunks[line[3]] = [line]
        elif line[3] == res_name:
            residue_chunks[line[3]].append(line)
        elif line[3] != res_name and res_name in residue_chunks.keys():
            break
        else:
            return {res_name: []}

    return residue_chunks


def discrimiate(residue, pdb_file):
    parsed, lines = read_pdb(f"./{pdb_file}")
    residue_chunk = get_residue_chunk(pdb_file, residue)
    atoms = [line[2] for line in residue_chunk[residue]]
    if "N6" in atoms:
        return f"The residue is a modification of ADENINE(A)"
    elif "O6" in atoms:
        return f"The residue is a modification of GUANINE(G)"
    elif "O4" in atoms:
        return f"The residue is a modification of URACIL(U)"
    elif "N4" in atoms:
        return f"The residue is a modification of CYTOCINE(C)"
    else:
        return "Couldn't figure it out :( "


print()
print("*********************************")
print(discrimiate("6M", "6ek0.pdb"))
print("*********************************")
print()
