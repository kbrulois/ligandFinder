



import numpy as np
import scipy.sparse
import sympy
import sklearn.datasets
import sklearn.feature_extraction.text
import umap
import umap.plot
import matplotlib.pyplot as plt
import sparse
import pandas as pd


dat = pd.read_csv('/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_lig1-15.csv')

dat = dat.iloc[:,77:]

dat = scipy.sparse.csr_matrix(dat)

umap_res = umap.UMAP(metric='euclidean', n_neighbors = 50, random_state=42, low_memory=True).fit(dat)

umap_embed = umap_res.embedding_

umap_embed = pd.DataFrame(umap_embed, columns=['UMAP1_eu_n30', 'UMAP2_eu_n30'])

umap_embed.to_csv('/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP_eu_n30.csv', index=False)




