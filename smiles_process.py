import lz4.frame
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from tucan.canonicalization import canonicalize_molecule
from tucan.io import graph_from_molfile_text
from tucan.serialization import serialize_molecule

fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
uncharger = rdMolStandardize.Uncharger()

def standardize(mol):
    clean_mol = rdMolStandardize.Cleanup(mol)
    bigger_clean = rdMolStandardize.FragmentParent(clean_mol)
    bigger_clean = uncharger.uncharge(bigger_clean)
    return bigger_clean

def fingerprint(mol):
    return fpgen.GetFingerprint(mol)

def process_smiles(inp):
    try:
        nid, smiles = inp
        mol = Chem.MolFromSmiles(smiles)
        smol = standardize(mol)
        if smol is not None:
            sim_fp = fingerprint(smol)
            sub_fp = Chem.PatternFingerprint(smol)
            return (nid, smiles, sim_fp, sub_fp)
        else:
            return None
    except:
        print("Failed to process", smiles)
        return None
