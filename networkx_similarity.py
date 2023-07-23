import pickle

import networkx as nx

from processor import microtucan, tucanize


def get_networkx_similarity(tucan: bytes) -> float:
    print(tucan.decode("utf8"))


#structure_db = pickle.load(open("data/structures_lz4.pkl", "rb"))

def get_bonds_from_tucan(tucan: str) -> list[tuple[int,int]]:
    bonds = tucan.split("/")[1]
    for bond in bonds.split(")"):
        atoms = bond[1:].split("-")
        if len(atoms) == 2:
            yield int(atoms[0]), int(atoms[1])

def get_graph_from_tucan(tucan: str) -> nx.Graph:
    g = nx.Graph()
    for bond in get_bonds_from_tucan(tucan):
        g.add_edge(bond[0], bond[1])
    return g

tucan=tucanize("CN(C)CCC1=CNC2=C1C(=CC=C2)").strip()

print(tucan)
print(len(tucan))

print(get_graph_from_tucan(tucan))

mt=microtucan(tucan)
print(len(mt))

