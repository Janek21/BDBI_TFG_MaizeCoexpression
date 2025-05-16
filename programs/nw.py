import pandas as pd

# Rows = genes, columns = samples
df = pd.read_csv("../data/wlen/data_wlen.csv", index_col=0)
##Filter by blue module

correlation_matrix = df.T.corr(method='pearson')


import numpy as np

threshold = 0.9  # You can adjust this
adjacency_matrix = (correlation_matrix.abs() > threshold).astype(int)
np.fill_diagonal(adjacency_matrix.values, 0)  # Remove self-loops


import networkx as nx

G = nx.from_pandas_adjacency(adjacency_matrix)


import matplotlib.pyplot as plt

plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, seed=42)
nx.draw_networkx(G, pos, node_size=50, with_labels=False, edge_color='gray')
plt.title("Gene Co-Expression Network")
plt.show()

