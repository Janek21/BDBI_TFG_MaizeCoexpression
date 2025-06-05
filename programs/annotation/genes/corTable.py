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


def main():
    filepath="./annotation/correlation_tables/brown_geneCorrelation.txt"
    #filepath="./annotation/correlation_tables/turquoise_geneCorrelation.txt"
    #Zm00001eb068950 #black
    #Zm00001eb255730 #brown
    #Zm00001eb000110 #darkred
    
    gint="Zm00001eb255730"
    corTable=pd.read_table(filepath, sep="\t")
    cor_edgelist=corrCleaner(corTable, 0.8)
    nw=nwLoad(cor_edgelist)
    nwShow(nw)#, "./graphFigs/fig1.png")

    print("Louvain:")
    louvain_comms = nx.community.louvain_communities(nw, weight="correlation", seed=42)
    
    print(f"Number of communities: {len(louvain_comms)}")

    # Draw each community
    for i, comm in enumerate(louvain_comms):
        if len(comm)>1:
            if gint in comm:
                print(comm)
                subgraph=nw.subgraph(comm)
                s_pos=nx.spring_layout(subgraph, seed=42)
                
                nx.draw(subgraph, pos=s_pos, with_labels=False, node_color=f"C{i}", node_size=100, font_size=6, edge_color="gray")
                
                #terms for the nodes and labels
                labels={node: node for node in subgraph.nodes if node==gint}
                nx.draw_networkx_labels(subgraph, pos=s_pos, labels=labels)
                nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=labels, node_color="black", node_size=100)
                
                plt.show()
                #plt.savefig(f"./graphFigs/{i}_fig.png")
                plt.close()
    


if __name__ == "__main__":
    main()