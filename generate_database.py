#!/usr/bin/env python3
import bz2
import multiprocessing
import pickle
import time
import csv

from concurrent.futures import ProcessPoolExecutor

from rdkit.Chem import AllChem

from smiles_process import process_smiles

start = time.time()

# read smiles.csv (two columns, structure_link and smiles)
with open("data/smiles.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    stuff = [x for x in reader]
    links, smileses = zip(*stuff)

max_workers = multiprocessing.cpu_count()


fps = []
processed_smileses = []
processed_links = []

with ProcessPoolExecutor(max_workers=max_workers) as executor:
    results = executor.map(process_smiles, enumerate(smileses), chunksize=16)
    for result in results:
        if result is not None:
            mid, smiles, fp = result
            fps.append(fp)
            processed_smileses.append(smiles)
            processed_links.append(links[mid].replace("http://www.wikidata.org/entity/", ""))


print(f"Precalculations time: {time.time() - start:.2f} s")
structure = {
    "smiles": processed_smileses,
    "links": processed_links,
    "fps": fps
}

with bz2.open("data/structure_fps.pkl.bz2", "wb") as f:
    f.write(pickle.dumps(structure))