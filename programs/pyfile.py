#!/usr/bin/env python3

print("in", flush=True)
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like

print("loaded", flush=True)
# Load Edge List CSV
df_edges = pd.read_csv("Cytoscape_edges.csv", sep='\t')
print(df_edges, flush=True)

# Load Node Colors
df_colors = pd.read_csv("Cytoscape_nodes.txt", sep='\t', header=0, names=['node', 'RM', 'color'])
print(df_colors, flush=True)
color_dict = dict(zip(df_colors["node"], df_colors["color"]))
print(df_colors.head(), flush=True)

# Create Graph
print("Creating graph", flush=True)
G = nx.from_pandas_edgelist(df_edges, source=df_edges.columns[0], target=df_edges.columns[1])

# Assign node colors
# Create a color dictionary
print("Creating color dictionary", flush=True)
color_dict = {node: color for node, color in zip(G.nodes(), df_colors['color'])}
# Assign colors to nodes, using a default color for any missing nodes
print("Assigning colors to nodes", flush=True)
node_colors = [color_dict.get(node, "gray") for node in G.nodes()]

print("Checking color validity...", flush=True)
valid_colors = []
for color in node_colors:
    if is_color_like(color):
        valid_colors.append(color)
    else:
        print(f"Warning: '{color}' is not a valid color. Using default color.", flush=True)
        valid_colors.append("gray")  # Use a default color for invalid entries

# Use valid_colors instead of node_colors in nx.draw()
print("drawing graph", flush=True)
nx.draw(G, with_labels=True, node_size=500, node_color=valid_colors, edge_color="gray", font_size=8)

# Save to PDF
print("Saving to pdf", flush=True)
plt.savefig("network.pdf", format="pdf")
plt.show()

print("Network saved as network.pdf", flush=True)

