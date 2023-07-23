import lz4.frame
from rdkit import Chem
from tucan.canonicalization import canonicalize_molecule
from tucan.io import graph_from_molfile_text
from tucan.serialization import serialize_molecule


def tucanize(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    m = Chem.MolToMolBlock(mol)
    molecule = graph_from_molfile_text(m)
    canonical_molecule = canonicalize_molecule(molecule)
    return serialize_molecule(canonical_molecule)


def microtucan(tucan: str) -> bytes:
    """Magic binary format"""
    out = bytes()
    split = tucan.split("/")
    out += split[0].encode("utf8") + b"/"
    bonds = split[1]
    out += "/".join(split[2:]).encode("utf8")
    for bond in bonds.split(")"):

        atoms = bond[1:].replace("(", "").split("-")
        if len(atoms) == 2:
            for atom in atoms:
                v = int(atom)
                if 255 < v < 65535:
                    out += v.to_bytes(3, "big")  # We prefix the longer ones with a 0 if the mol has more than 255 atoms
                else:
                    out += v.to_bytes(1)

    return out


def process_smiles(inp):
    try:
        nid, smiles = inp
        tucan_string = microtucan(tucanize(smiles))
        return (nid, smiles, tucan_string,
                len(lz4.frame.compress(tucan_string)))
    except:
        print("Failed to process", smiles)
        return None
