#!/usr/bin/env python3
import multiprocessing
import pickle
import time
import csv

from concurrent.futures import ProcessPoolExecutor

from smiles_process import process_smiles

start = time.time()

# read smiles.csv (two columns, structure_link and smiles)
with open("data/smiles.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    stuff = [x for x in reader]
    links, smileses = zip(*stuff)

max_workers = multiprocessing.cpu_count()
tucans = []
tucanlz4s = []
processed_smileses = []
processed_links = []

with ProcessPoolExecutor(max_workers=max_workers) as executor:
    results = executor.map(process_smiles, enumerate(smileses), chunksize=16)
    for result in results:
        if result is not None:
            mid, smiles, tucan, tucanlz4 = result
            tucans.append(tucan)
            tucanlz4s.append(tucanlz4)
            processed_smileses.append(smiles)
            processed_links.append(links[mid].replace("http://www.wikidata.org/entity/", ""))


print(f"Precalculations time: {time.time() - start:.2f} s")
structure = {
    "smiles": processed_smileses,
    "links": processed_links,
    "tucans": tucans,
    "tucanlz4s": tucanlz4s
}

pickle.dump(structure, open("data/structures_lz4.pkl", "wb"))
