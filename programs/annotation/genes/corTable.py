#!/usr/bin/env python3

#Load gene Correlation table ino networkX, genes are nodes, correlations above x value are edges

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

def corrCleaner(table, ths=0.7):
    #Join the table as gene to gene relations: rows are [gene1 gene2 correlationCoeff]
    print(type(table))
    corList=table.stack().reset_index()
    corList.columns=["G1", "G2", "coefficient"]
    print(len(corList))
    
    #Take out the self-loops
    corList_noSL=corList[corList["G1"]!=corList["G2"]]
    print(len(corList_noSL))
    #filter by threshold (observed in genCor_clustering script(plot))
    correlations=corList_noSL[corList_noSL["coefficient"].abs()>=ths]
    print(len(correlations))
    
    return correlations

def nwLoad(edgelist):
    nw=nx.from_pandas_edgelist(edgelist, "G1", "G2", edge_attr="coefficient")
    return nw

def nwShow(nw, filename=None):
    nx.draw(nw,node_size=10)
    
    #Show or save
    if filename!=None:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close()
    return

def main():
    filepath="./annotation/ME_correlation/lightcyan_geneCorrelation.txt"
    corTable=pd.read_table(filepath, sep="\t")
    cor_edgelist=corrCleaner(corTable, 0.6)
    nw=nwLoad(cor_edgelist)
    nwShow(nw, "./graphFigs/fig1.png")

    print("Louvain:")
    louvain_comms = nx.community.louvain_communities(nw, weight="correlation", seed=42)
    for x in louvain_comms:
        print(len(x))

    # Draw each community
    for i, comm in enumerate(louvain_comms):
        if len(comm)>1:
            subgraph=nw.subgraph(comm)
            
            nx.draw(subgraph,with_labels=True, node_color=f"C{i}", node_size=100, font_size=6, edge_color="gray")

            plt.savefig(f"./graphFigs/{i}_fig.png")
            plt.close()
    
        

if __name__ == "__main__":
    main()