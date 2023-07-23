 This application is available on [Streamlit](https://lotus-tanimoto.streamlit.app/)

I had a really quick and dirty experiment for using lz4 and Tucan to find similar molecules, this is the version
using standard tanimoto and fingerprints from Rdkit.

The code of this is on my alter-github account: https://github.com/bjonnh-streamlit/mol_search

The dataset is from the [LOTUS](https://lotus.nprod.net/) database and [Wikidata](https://www.wikidata.org).

The idea to use Tucan and play with all that came from [Adriano Rutz](https://adafede.github.io/), we are both
part of the team behind LOTUS.

It is using [Streamlit](https://streamlit.io)  for its web ui.

[Rdkit](https://www.rdkit.org) is used for the molecule massaging.

See those other projects for more information on why I am trying that (and Daniel's work is what started me on that
contraption):

- https://github.com/daenuprobst/molzip
- https://arxiv.org/abs/2212.09410


To reproduce:
- Install and run poetry update
- Run `get_smiles_wikidata.py`  (takes 12s)
- Run `generate_database.py`    (takes 1min on a Ryzen 7 3700X)
- Run `streamlit run main.py`   (almost instant)


## **Me**

- Personal website: https://www.bjonnh.net
- You can find me on mastodon too(t): https://mastodon.social/@bjonnh


## **Data safety**

Your molecules are never stored, the only people that could eventually see them are the streamlit people. 
And I doubt they care about your molecules. But if your molecules are super secret, it is like your extremities, don't insert 
them in machines you don't know or understand.

Also, this is an experimental tool meant to test things, you're not supposed to rely on it for anything important, and
it comes with no warranty or support whatsoever.

## **License and legalese**

https://raw.githubusercontent.com/bjonnh-streamlit/mol_search/main/LICENSE
