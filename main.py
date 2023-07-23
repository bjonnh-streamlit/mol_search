#!/usr/bin/env python3
import base64
import pickle
import time

import lz4.frame
import numpy as np
import streamlit as st
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from tucan.canonicalization import canonicalize_molecule
from tucan.io import graph_from_molfile_text
from tucan.serialization import serialize_molecule

start = time.time()

structure_db = pickle.load(open("data/structures_lz4.pkl", "rb"))
smiles = structure_db["smiles"]
tucans = structure_db["tucans"]
tucanlz4s = structure_db["tucanlz4s"]
links = structure_db["links"]

end_dl_time = time.time()


def search(smiles):
    mol = Chem.MolFromSmiles(smiles)
    m = Chem.MolToMolBlock(mol)
    molecule = graph_from_molfile_text(m)
    canonical_molecule = canonicalize_molecule(molecule)
    tucan1 = serialize_molecule(canonical_molecule).encode("utf-8")
    tucanlz4_1 = len(lz4.frame.compress(tucan1))

    out = []
    for j in range(len(tucans)):
        tucan2 = tucans[j]
        tucanlz4_2 = tucanlz4s[j]

        lc = len(lz4.frame.compress((tucan1 + tucan2)))
        tucanlz4_score = (lc - min(tucanlz4_2, tucanlz4_1)) / max(tucanlz4_2, tucanlz4_1)

        out.append((j, tucanlz4_score))
    return out


def molecule_svg(smiles):
    d2d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    d2d.DrawMolecule(Chem.MolFromSmiles(smiles))
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def render_svg(svg):
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

st.title("LOTUS LZ4 searcher")
with st.expander("About"):
    st.markdown(open("README.md").read())

query = st.text_input(label="SMILES (short ones work really really badly you've been warned)", value="O=C1c3c(O/C(=C1/O)c2ccc(O)c(O)c2)cc(O)cc3O")
level = st.slider("Level (lower means tighter search)", min_value=0.0, max_value=0.8, value=0.3, step=0.01)
if query != "":
    start = time.time()
    try:
        render_svg(molecule_svg(query))
        scores = [score for score in search(query)]

        # with transparency
        js, tucanlz4_scores = zip(*scores)

        count = 0
        results = []
        for idx, v in enumerate(tucanlz4_scores):
            if v < level:
                results.append([v, idx])
        sorted_results = sorted(results, key=lambda x: x[0])
        st.header("Results")

        arr = np.array(tucanlz4_scores)
        fig, ax = plt.subplots()
        ax.hist(arr, bins=20)
        st.write("Histogram of LZ4 scores (you want to put your level way before the peak on the left!")
        st.pyplot(fig)
        for result in sorted_results:
            st.divider()
            m = smiles[result[1]]
            render_svg(molecule_svg(m))
            st.markdown(f"[Link to Wikidata]({links[result[1]]})")
            st.write(f"{m} - score: {result[0]}")

        st.write(
            f"Search elapsed time: {time.time() - start:.2f}s (data_loading: {end_dl_time - start:.2f}s) with {count} matches < {level}")
    except:
        st.error("Your molecule is likely invalid.")