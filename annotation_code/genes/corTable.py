#!/usr/bin/env python3

#Load gene Correlation table ino networkX, genes are nodes, correlations above x value are edges

import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use('gtk3agg') 
import matplotlib.pyplot as plt

def corrCleaner(table, ths=0.7):
    #Join the table as gene to gene relations: rows are [gene1 gene2 correlationCoeff]
    corList=table.stack().reset_index()
    corList.columns=["G1", "G2", "coefficient"]
    
    #Take out the self-loops
    corList_noSL=corList[corList["G1"]!=corList["G2"]]
    #filter by threshold (observed in genCor_clustering script(plot))
    correlations=corList_noSL[corList_noSL["coefficient"].abs()>=ths]
    
    return correlations

def nwLoad(edgelist):
    nw=nx.from_pandas_edgelist(edgelist, "G1", "G2", edge_attr="coefficient")
    return nw

def nwShow(nw, gr_title="No module", filename=None):
    plt.figure()  
    plt.title(gr_title)
    nx.draw(nw,node_size=10)
    #Show or save
    if filename!=None:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close()

def subgraph_show(subgraph, comm_index, gene_of_interest="", folderpath=None):
    s_pos=nx.spring_layout(subgraph, seed=42)
    
    nx.draw(subgraph, pos=s_pos, with_labels=False, node_color=f"C{comm_index}", node_size=100, font_size=6, edge_color="gray")
    
    #terms for the nodes and labels
    labels={node: node for node in subgraph.nodes if node==gene_of_interest}
    nx.draw_networkx_labels(subgraph, pos=s_pos, labels=labels)
    nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=labels, node_color="black", node_size=100)
    
    if folderpath!=None:
        plt.savefig(f"{folderpath}/{comm_index}_fig.png")
    else:
        plt.show()
    
    
import GO_comm as gc
from collections import defaultdict

def main():
    filepath="./correlation_tables/brown_geneCorrelation.txt"
    #filepath="./correlation_tables/turquoise_geneCorrelation.txt"
    #Zm00001eb068950 #black
    #Zm00001eb255730 #brown
    #Zm00001eb000110 #darkred
    
    gint="Zm00001eb255730"
    gint=""
    corTable=pd.read_table(filepath, sep="\t")
    cor_edgelist=corrCleaner(corTable, 0.6)
    #create the graph
    nw=nwLoad(cor_edgelist)
    nwShow(nw)#, "./graphFigs/fig1.png")

    #search communities
    print("Louvain:")
    louvain_comms = nx.community.louvain_communities(nw, weight="correlation", seed=42)
    
    print(f"Number of communities: {len(louvain_comms)}")

    #draw each community
    for i, comm in enumerate(louvain_comms):
        if len(comm)>1:
            #if gint in comm:
                comm=list(comm)
                subgraph=nw.subgraph(comm)
                
                subgraph_show(subgraph, i, gene_of_interest=gint) #, folderpath="./graphFigs")
                #plt.close()

                annot=gc.functional_analysis(comm, "zmays")

                gc.genesToGo(len(comm), len(annot))
                plt.close()
## return annotated genes and color them in plot


if __name__ == "__main__":
    main()