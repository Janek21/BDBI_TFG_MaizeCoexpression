#!/usr/bin/env python3

#Load gene Correlation table ino networkX, genes are nodes, correlations above x value are edges

import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use('gtk3agg') 
import matplotlib.pyplot as plt

import GOannotation as gan

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

def nwShow(nw, dimlist, gr_title, filename=None):
    plt.figure(figsize=(dimlist[0], dimlist[1]), dpi=dimlist[2])
    plt.title(gr_title)
    nx.draw(nw, node_size=10, node_color=gr_title)
    #Show or save
    if filename!=None:
        plt.savefig(filename, dpi=dimlist[2], bbox_inches="tight")
    else:
        plt.show()
    plt.close()
    
def tableWriter(i, functions_dic, openfile):
    ukg_list=functions_dic.keys()
    
    for unknownGene in ukg_list: #1 row per unknown gene in the comm
        kg=functions_dic[unknownGene][0] #list of known genes
        funct=functions_dic[unknownGene][1] #list of functions
        print(f"C{i}\t{unknownGene}\t{'.'.join(kg)}\t{'.'.join(funct)}", file=openfile)
        
    return


def main():
    module="turquoise"
    filepath=f"./correlation_tables/{module}_geneCorrelation.txt"
    
    #create network
    corTable=pd.read_table(filepath, sep="\t")
    cor_edgelist=corrCleaner(corTable, 0.6)
    nw=nwLoad(cor_edgelist)
    
    #set plot size
    width_px=1920
    height_px=1120
    dpi=100
    width_in=width_px / dpi
    height_in=height_px / dpi
    dimlist=[width_in, height_in, dpi]
    nwShow(nw, dimlist, module, f"./graphFigs/{module}_general.png") #draw network
    
    #save gene information
    commTablepath=f"./genes/community_tables/{module}_commF.txt"
    comT=open(commTablepath, "w")
    print(f"Community\tUnknown\Annotated\tFunctions", file=comT) #header
    
    
    #Communities
    louvain_comms=nx.community.louvain_communities(nw, weight="correlation", seed=42)
    print(f"Number of communities: {len(louvain_comms)}")
    
    #
    color_map={}
    colors=plt.cm.get_cmap('viridis', len(louvain_comms))
    
    #Draw each community
    for i, comm in enumerate(louvain_comms):
        if len(comm)>2:
            #if i==5: #look into a particular community
                fd=gan.main(module, comm, 10, 10) #functions dictionary
                subgraph=nw.subgraph(comm)
                s_pos=nx.spring_layout(subgraph, seed=42)
                
                plt.figure(figsize=(width_in, height_in), dpi=dpi)
                plt.title(f"Community {i}: {len(comm)} genes\n Known: {len(comm)-len(fd.keys())} Unknown:{len(fd.keys())}")
                nx.draw(subgraph, pos=s_pos, with_labels=False, node_size=100, edge_color="gray", node_color=module)
                nx.draw_networkx_nodes(subgraph, pos=s_pos, node_size=100, linewidths=1, edgecolors="black")
                
                #terms for the nodes and labels
                nestList_cg=[val[0] for val in fd.values()] #get only corr genes
                corrGen=sum(nestList_cg, []) #unnest the list
                
                unknown_genes={node: node for node in subgraph.nodes if node in fd.keys()}
                #k_genes={node: node for node in subgraph.nodes if node in corrGen} #find out closely correlated genes (find out if a gene in a community has no close correls)
                nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=unknown_genes, node_color="black", node_size=100)
                #nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=k_genes, node_color=f"red", node_size=100) #view closerly correlated genes
                
                plt.savefig(f"./graphFigs/{module}_{i}.png", dpi=dpi, bbox_inches="tight") #save plot
                #plt.show()
                plt.close()
                
                #save info to table for easy parsing
                tableWriter(i, fd, comT)
                
                #grep C1 turquoise_commF.txt |cut -f4 -d$'\t'|cut -d'.' -f1|sort|uniq -c|sort -n
                
                
                #classify nodes by community
        for node in comm:
            print(node)
            color_map[node]=colors(i)
                
    comT.close()
    
    node_col=[color_map[node] for node in nw.nodes()]
    plt.figure(figsize=(dimlist[0], dimlist[1]), dpi=dimlist[2])
    nw_pos=nx.spring_layout(nw, seed=42)
    nx.draw(nw,node_size=10)
    nx.draw_networkx_nodes(nw, pos=nw_pos, node_colors=node_col, node_size=100, linewidths=1, edgecolors="black")
    plt.show()
    


if __name__ == "__main__":
    main()