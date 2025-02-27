import pandas as pd
import numpy as np
import networkx as nx

# Load the gene expression data
data = pd.read_csv("/home/janek/Desktop/Intership/work/data/DK105.csv", index_col=0)

print("e")
# Transpose the data if genes are in columns
data = data.T

# Compute correlation matrix
cor_matrix = data.corr(method='pearson')

# Set threshold for network edges
threshold = 0.7
edges = np.where(abs(cor_matrix) > threshold)

# Create graph
G = nx.Graph()
for i, j in zip(edges[0], edges[1]):
    if i != j:  # Avoid self-loops
        G.add_edge(data.columns[i], data.columns[j], weight=cor_matrix.iloc[i, j])

# Save network
nx.write_graphml(G, "gene_network.graphml")

# Optional: Draw the network
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 10))
nx.draw(G, node_size=5, edge_color="gray", with_labels=False)
plt.show()
