#!/usr/bin/env python3
import base64
import pickle
import time

import streamlit as st
from matplotlib import pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary
from rdkit.Chem.Draw import rdMolDraw2D

from smiles_process import fingerprint, standardize

start = time.time()
st.set_page_config(page_title="LOTUS searcher", page_icon=":lotus:", layout="centered",
                   initial_sidebar_state="auto", menu_items=None)


@st.cache_resource(ttl=3600)
def load_data():
    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    fps = rdSubstructLibrary.PatternHolder()
    with open("data/structure_fps.pkl","rb") as f:
        data = pickle.load(f)
    for idx in range(len(data["sub_fps"])):
        mols.AddSmiles(data["smileses"][idx])
        fps.AddFingerprint(data["sub_fps"][idx])

    data["library"] = rdSubstructLibrary.SubstructLibrary(mols, fps)
    return data


@st.cache_data(ttl=3600)
def readme():
    return "".join([line for line in open("README.md").readlines() if not line.startswith(" ")])


structure_db = load_data()


def search(fp: bytes) -> list[tuple[int, float]]:
    results = DataStructs.BulkTanimotoSimilarity(fp, structure_db["sim_fps"])
    return [(j, score) for j, score in enumerate(results)]

def ss_search(fp, mol) -> list[tuple[int, float]]:
    out = []
    for j in structure_db["library"].GetMatches(mol):
        out.append((j, DataStructs.TanimotoSimilarity(fp, structure_db["sim_fps"][j])))
    return out



def molecule_svg(mol):
    d2d = rdMolDraw2D.MolDraw2DSVG(250, 200)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()


def render_svg(svg):
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)


st.title("LOTUS searcher")
with st.expander("About"):
    st.markdown(readme())

st.write("Choose an example, or type your SMILES below")
c1, c2, c3, c4 = st.columns(4)

amarogentin = "C=C[C@@H]1[C@@H]2CCOC(=O)C2=CO[C@H]1O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)OC(=O)C4=C(C=C(C=C4C5=CC(=CC=C5)O)O)O"

if "input_query" not in st.session_state:
    st.session_state["input_query"] = amarogentin

if c1.button("Amarogentin"):
    st.session_state["input_query"] = amarogentin
    st.experimental_rerun()
if c2.button("Quassin"):
    st.session_state[
        "input_query"] = "C[C@@H]1C=C(C(=O)[C@]2([C@H]1C[C@@H]3[C@@]4([C@@H]2C(=O)C(=C([C@@H]4CC(=O)O3)C)OC)C)C)OC"
    st.experimental_rerun()
if c3.button("Absinthin"):
    st.session_state[
        "input_query"] = "C[C@H]1[C@@H]2CC[C@]([C@@H]3[C@H]4[C@H]5C=C([C@@]6([C@H]4C(=C3[C@H]2OC1=O)C)[C@@H]5[C@@](CC[C@@H]7[C@@H]6OC(=O)[C@H]7C)(C)O)C)(C)O"
    st.experimental_rerun()
if c4.button("Quinine"):
    st.session_state["input_query"] = "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O"
    st.experimental_rerun()

query = st.text_input(label="SMILES",
                      key="input_query")
ss_mode = st.checkbox("Sub-structure search")
if query:
    try:
        mol = standardize(Chem.MolFromSmiles(query))
        render_svg(molecule_svg(mol))
        fp = fingerprint(mol)

        if ss_mode:
            scores = ss_search(fp, mol)
        else:
            scores = search(fp)
            fig, ax = plt.subplots()
            ax.hist([score[1] for score in scores], bins=20)
            st.write("Histogram of tanimoto scores (you want to put your level way before the peak on the left!")

            st.pyplot(fig, use_container_width=True)
            level = st.slider("Level (higher means tighter search)", min_value=0.0, max_value=1.0, value=0.8, step=0.01)
            scores = [score for score in scores if score[1] >= level]
        st.header("Results")
        st.write(
            "The free plan of streamlit is a bit memory limited (but it is still the best out there with its 1GB)"
            " so we will only show you 100 matches. If you run it locally you can see all of them.")

        start = time.time()

        count = len(scores)
        scores = sorted(scores, reverse=True, key=lambda x: x[1])

        scores = scores[0:100]
        st.write(f"Search time: {time.time() - start:.2f}s with {count} matches")

        for result in scores:
            st.divider()
            m = structure_db["smileses"][result[0]]
            render_svg(molecule_svg(structure_db["library"].GetMol(result[0])))
            wid = structure_db["links"][result[0]]
            st.markdown(f"[Link to Wikidata](http://www.wikidata.org/entity/{wid})")

            st.text(m)

            st.progress(result[1], text="Tanimoto similarity: {:.2f}".format(result[1]))

    except Exception as e:
        st.error(f"Your molecule is likely invalid. {e}")
