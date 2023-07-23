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

st.set_page_config(page_title="LOTUS LZ4 searcher", page_icon=":lotus:", layout="centered",
                   initial_sidebar_state="auto", menu_items=None)


@st.cache_data(ttl=3600)
def load_data():
    return pickle.load(open("data/structures_lz4.pkl", "rb"))

@st.cache_data(ttl=3600)
def readme():
    return "".join([line for line in open("README.md").readlines() if not line.startswith(" ")])

structure_db = load_data()


def search(smiles):
    mol = Chem.MolFromSmiles(smiles)
    m = Chem.MolToMolBlock(mol)
    molecule = graph_from_molfile_text(m)
    canonical_molecule = canonicalize_molecule(molecule)
    tucan1 = serialize_molecule(canonical_molecule).encode("utf-8")
    tucanlz4_1 = len(lz4.frame.compress(tucan1))

    out = []
    for j in range(len(structure_db["tucans"])):
        tucan2 = structure_db["tucans"][j]
        tucanlz4_2 = structure_db["tucanlz4s"][j]

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
    st.markdown(readme())


st.write("Choose an example, or type your SMILES below")
c1, c2, c3, c4 = st.columns(4)

amarogentin = "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O"

if "input_query" not in st.session_state:
    st.session_state["input_query"] = amarogentin

if c1.button("Amarogentin"):
    st.session_state["input_query"] = amarogentin
if c2.button("Quassin"):
    st.session_state["input_query"] = "C[C@@H]1C=C(C(=O)[C@]2([C@H]1C[C@@H]3[C@@]4([C@@H]2C(=O)C(=C([C@@H]4CC(=O)O3)C)OC)C)C)OC"
if c3.button("Absinthin"):
    st.session_state["input_query"] = "C[C@H]1[C@@H]2CC[C@]([C@@H]3[C@H]4[C@H]5C=C([C@@]6([C@H]4C(=C3[C@H]2OC1=O)C)[C@@H]5[C@@](CC[C@@H]7[C@@H]6OC(=O)[C@H]7C)(C)O)C)(C)O"
if c4.button("Quinine"):
    st.session_state["input_query"] = "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O"

query = st.text_input(label="SMILES (short ones work really really badly you've been warned)",
                            key="input_query")



try:
    render_svg(molecule_svg(query))
    scores = search(query)
    fig, ax = plt.subplots()
    ax.hist([score[1] for score in scores], bins=20)
    st.write("Histogram of LZ4 scores (you want to put your level way before the peak on the left!")

    st.pyplot(fig, use_container_width=True)
    level = st.slider("Level (lower means tighter search)", min_value=0.0, max_value=0.8, value=0.3, step=0.01)

    st.write(
        "The free plan of streamlit is a bit memory limited (but it is still the best out there with its 1GB)"
        " so we will only show you 100 matches. If you run it locally you can see all of them.")

    start = time.time()


    st.header("Results")
    results = [score for score in scores if score[1] < level]
    count = len(results)
    sorted_results = sorted(results, key=lambda x: x[1])[0:100]
    st.write(
        f"Search time: {time.time() - start:.2f}s with {count} matches")

    for result in sorted_results:
        st.divider()
        m = structure_db["smiles"][result[0]]
        render_svg(molecule_svg(m))
        wid = structure_db["links"][result[0]]
        st.markdown(f"[Link to Wikidata](http://www.wikidata.org/entity/{wid})")

        st.text(m)

        st.progress(1 - result[1], text="Distance: {:.2f}".format(result[1]))

except Exception as e:
    st.error(f"Your molecule is likely invalid.")
