import gzip
import struct

import lz4.frame
from rdkit import Chem, DataStructs
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

def compute_scores(indices):
    block_start, size, data_block = indices
    tucangzip = data_block["tucangzips"]
    tucans = data_block["tucans"]
    rdkitfps = data_block["rdkitfps"]

    blk_total = len(tucans)
    content = bytes()
    for i in range(block_start, min(block_start + size, blk_total)):
        tucan1 = tucans[i]
        tucangzip1 = tucangzip[i]
        rdkitfp1 = rdkitfps[i]
        for j in range(i + 1, blk_total):
            tucan2 = tucans[j]
            tucangzip2 = tucangzip[j]
            rdkitfp2 = rdkitfps[j]

            tanimoto_score = DataStructs.FingerprintSimilarity(rdkitfp1, rdkitfp2)
            lc = len(gzip.compress((tucan1 + tucan2)))
            tucangzip_score = (lc - min(tucangzip2, tucangzip1)) / max(tucangzip2, tucangzip1)
            content += (i.to_bytes(4) + j.to_bytes(4) + struct.pack("!f", tanimoto_score) +
                        struct.pack("!f", tucangzip_score))

    with open(f"data/scores_{block_start}_l_{size}.bin", "wb") as f:
        f.write(content)



def compute_scores_search(indices):
    start, blks, data_block = indices
    tucans = data_block["tucans"]
    tucangzips = data_block["tucangzips"]
    tucanlz4s = data_block["tucanlz4s"]
    rdkitfps = data_block["rdkitfps"]
    tucan1 = data_block["tucan1"]
    tucangzip1 = data_block["tucangzip1"]
    tucanlz4_1 = data_block["tucanlz4_1"]
    rdkitfp1 = data_block["rdkitfp1"]
    out = []
    for j in range(len(tucans)):
        tucan2 = tucans[j]
        tucangzip2 = tucangzips[j]
        tucanlz4_2 = tucanlz4s[j]
        rdkitfp2 = rdkitfps[j]

        tanimoto_score = DataStructs.FingerprintSimilarity(rdkitfp1, rdkitfp2)

        lc = len(gzip.compress((tucan1 + tucan2)))
        tucangzip_score = (lc - min(tucangzip2, tucangzip1)) / max(tucangzip2, tucangzip1)

        lc = len(lz4.frame.compress((tucan1 + tucan2)))
        tucanlz4_score = (lc - min(tucanlz4_2, tucanlz4_1)) / max(tucanlz4_2, tucanlz4_1)

        out.append((j, tanimoto_score, tucangzip_score, tucanlz4_score))
    return out
