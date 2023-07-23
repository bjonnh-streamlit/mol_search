#!/usr/bin/env python3
import gzip
import pickle

import lz4.frame

### Searching for molecules that can't find themselves

structure_db = pickle.load(open("data/structures_lz4.pkl", "rb"))


def lz4_compress(data: bytes, level: int = 3) -> bytes:
    return lz4.frame.compress(data, compression_level=level)


def get_tucan_lz4_score(tucan: bytes) -> float:
    tucanlz4 = len(lz4_compress(tucan))
    lc = len(lz4_compress((tucan + tucan)))
    return lc / tucanlz4 - 1

def get_tucan_highlz4_score(tucan: bytes) -> float:
    tucanlz4 = len(lz4_compress(tucan,16))
    lc = len(lz4_compress(tucan + tucan,16))
    return lc / tucanlz4 - 1


def get_tucan_gzip_score(tucan: bytes) -> float:
    tucanlz4 = len(gzip.compress(tucan))
    lc = len(gzip.compress((tucan + tucan)))
    return lc / tucanlz4 - 1


for idx, tucan in enumerate(structure_db["tucans"]):
    tucanlz4_score = get_tucan_lz4_score(tucan)
    tucanhighlz4_score = get_tucan_highlz4_score(tucan)
    tucangzip_score = get_tucan_gzip_score(tucan)
    #networkx_similarity = get_networkx_similarity(tucan)

    if tucanlz4_score > 0.2:
        print(
            f"self_tucan_lz4={tucanlz4_score:.3f} self_tucan_highlz4={tucanhighlz4_score:.3f} self_tucan_gzip={tucangzip_score:.3f} {structure_db['smiles'][idx]}")
