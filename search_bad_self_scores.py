#!/usr/bin/env python3
import pickle

import lz4.frame

### Searching for molecules that can't find themselves

structure_db = pickle.load(open("data/structures_lz4.pkl", "rb"))

for idx, tucan in enumerate(structure_db["tucans"]):
    tucanlz4 = len(lz4.frame.compress(tucan))

    lc = len(lz4.frame.compress((tucan + tucan)))
    tucanlz4_score = lc/tucanlz4 - 1

    if tucanlz4_score > 0.25:
        print(f"{tucanlz4_score:.3f}: {structure_db['smiles'][idx]}")


