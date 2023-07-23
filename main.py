#!/usr/bin/env python3
import base64
import pickle
import time

import lz4.frame
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


st.set_page_config(page_title="LOTUS LZ4 searcher", page_icon=":lotus:", layout="centered",
                   initial_sidebar_state="auto", menu_items=None)

st.title("LOTUS LZ4 searcher")
with st.expander("About"):
    st.markdown(open("README.md").read())

query = st.text_input(label="SMILES (short ones work really really badly you've been warned)",
                      value="O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12")
level = st.slider("Level (lower means tighter search)", min_value=0.0, max_value=0.8, value=0.3, step=0.01)
st.write(
    "The free plan of streamlit is a bit memory limited so we will only show you 100 matches. If you run it locally you can see all of them.")
if query != "":
    start = time.time()
    try:
        render_svg(molecule_svg(query))
        scores = search(query)

        results = [score for score in scores if score[1] < level]
        count = len(results)
        sorted_results = sorted(results, key=lambda x: x[1])[0:100]
        st.header("Results")
        st.write(
            f"Search time: {time.time() - start:.2f}s with {count} matches")
        fig, ax = plt.subplots()
        ax.hist([score[1] for score in scores], bins=20)
        st.write("Histogram of LZ4 scores (you want to put your level way before the peak on the left!")

        st.pyplot(fig, use_container_width=True)

        for result in sorted_results:
            st.divider()
            m = smiles[result[0]]
            render_svg(molecule_svg(m))
            st.markdown(f"[Link to Wikidata](http://www.wikidata.org/entity/{links[result[0]]})")

            st.text(m)

            st.progress(1 - result[1], text="Distance: {:.2f}".format(result[1]))

    except Exception as e:
        st.error(f"Your molecule is likely invalid.")
