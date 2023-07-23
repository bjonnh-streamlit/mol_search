#!/usr/bin/env python3
import lz4.frame
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
    links, isomeric_smileses, canonical_smileses = zip(*stuff)
    smileses = []
    for i in range(len(isomeric_smileses)):
        if isomeric_smileses[i] == "":
            smileses.append(canonical_smileses[i])
        else:
            smileses.append(isomeric_smileses[i])

max_workers = multiprocessing.cpu_count()


p_smileses = []
p_sim_fps = []
p_sub_fps = []
p_links = []

with ProcessPoolExecutor(max_workers=max_workers) as executor:
    results = executor.map(process_smiles, enumerate(smileses), chunksize=1000)
    for result in results:
        if result is not None:
            mid, smiles, sim_fp, sub_fp = result
            p_smileses.append(smiles)
            p_sim_fps.append(sim_fp)
            p_sub_fps.append(sub_fp)
            p_links.append(links[mid].replace("http://www.wikidata.org/entity/", ""))


structure = {
    "smileses": p_smileses,
    "links": p_links,
    "sim_fps": p_sim_fps,
    "sub_fps": p_sub_fps
}


pickle.dump(structure, open("data/structure_fps.pkl", "wb"))
print(f"Precalculations time: {time.time() - start:.2f} s")
