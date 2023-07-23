import lz4.frame
from rdkit import Chem
from rdkit.Chem import AllChem
from tucan.canonicalization import canonicalize_molecule
from tucan.io import graph_from_molfile_text
from tucan.serialization import serialize_molecule

fpgen = AllChem.GetRDKitFPGenerator()
def process_smiles(inp):
    try:
        nid, smiles = inp
        mol = Chem.MolFromSmiles(smiles)
        fp = fpgen.GetFingerprint(mol)
        return (nid, smiles, fp)
    except:
        print("Failed to process", smiles)
        return None
