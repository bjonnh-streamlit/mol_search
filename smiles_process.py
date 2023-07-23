import lz4.frame
from rdkit import Chem
from tucan.canonicalization import canonicalize_molecule
from tucan.io import graph_from_molfile_text
from tucan.serialization import serialize_molecule


def process_smiles(inp):
    try:
        nid, smiles = inp
        mol = Chem.MolFromSmiles(smiles)
        m = Chem.MolToMolBlock(mol)
        molecule = graph_from_molfile_text(m)
        canonical_molecule = canonicalize_molecule(molecule)
        tucan_string = serialize_molecule(canonical_molecule).encode("utf-8")
        return (nid, smiles, tucan_string,
                len(lz4.frame.compress(tucan_string)))
    except:
        print("Failed to process", smiles)
        return None
