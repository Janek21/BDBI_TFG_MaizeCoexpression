import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load Edge List CSV
df_edges = pd.read_csv("your_file.csv")

# Load Node Colors
df_colors = pd.read_csv("cy_nodes.txt", header=None, names=["node", "color"])
color_dict = dict(zip(df_colors["node"], df_colors["color"]))

# Create Graph
G = nx.from_pandas_edgelist(df_edges, source=df_edges.columns[0], target=df_edges.columns[1])

# Assign node colors
node_colors = [color_dict.get(node, "gray") for node in G.nodes()]

# Draw the graph
plt.figure(figsize=(10, 10))
nx.draw(G, with_labels=True, node_size=500, node_color=node_colors, edge_color="gray", font_size=8)

# Save to PDF
plt.savefig("network.pdf", format="pdf")
plt.show()

print("Network saved as network.pdf")

