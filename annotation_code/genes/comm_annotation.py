#!/usr/bin/env python3

#Load gene Correlation table ino networkX, create Louvain communities and annotate the unknwn genes in this communities
#create tables with all genes annotated(known and unknown), for each inputed module

import matplotlib.patches
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use('gtk3agg') 
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like

import GOannotation as gan

def corrCleaner(table, ths=0.7):
    #Join the table as gene to gene relations: rows are [gene1 gene2 correlationCoeff]
    corList=table.stack().reset_index()
    corList.columns=["G1", "G2", "coefficient"]
    
    #Take out the self-loops
    corList_noSL=corList[corList["G1"]!=corList["G2"]]
    #filter by threshold (observed in geneClustering.Rmd)
    correlations=corList_noSL[corList_noSL["coefficient"].abs()>=ths]
    
    return correlations

def nwLoad(edgelist):
    nw=nx.from_pandas_edgelist(edgelist, "G1", "G2", edge_attr="coefficient")
    return nw

def nwShow(nw, dimlist, gr_title, filename=None):
    plt.figure(figsize=(dimlist[0], dimlist[1]), dpi=dimlist[2])
    plt.title(gr_title)
    if is_color_like(gr_title):
        nx.draw(nw, node_size=10, node_color=gr_title)
    else:
        nx.draw(nw, node_size=10)
    #Show or save
    if filename!=None:
        plt.savefig(filename, dpi=dimlist[2], bbox_inches="tight")
    else:
        plt.show()
    plt.close()
    
def tableWriter(i, functions_dic, openfile):
    gene_list=functions_dic.keys() #all genes
    known_G=[] #known gene list
    unk_G=[] #unknown gene list
    
    for g in gene_list: #1 row per unknown gene in the comm
        kg=functions_dic[g][0] #list of correlated(known) genes (or "Known" when gene is known)
        funct=functions_dic[g][1] #list of functions
        print(f"C{i}\t{g}\t{'.'.join(kg)}\t{'.'.join(funct)}", file=openfile)

        if "Known" in kg: #known gene list
            known_G.append(g)
        else: #unknown gene list
            unk_G.append(g)

    return [unk_G, known_G] #Unknown, known genes


def main(module, thrs=0.6, community_limit=2):
    filepath=f"./correlation_tables/{module}_geneCorrelation.txt"
    
    #create network
    corTable=pd.read_table(filepath, sep="\t")
    cor_edgelist=corrCleaner(corTable, thrs)
    nw=nwLoad(cor_edgelist)
    
    #set plot size
    width_px=1920
    height_px=1120
    dpi=100
    width_in=width_px / dpi
    height_in=height_px / dpi
    dimlist=[width_in, height_in, dpi]
    nwShow(nw, dimlist, module, f"./genes/graphFigs/{module}_general.png") #draw network
    
    #save gene information
    commTablepath=f"./genes/community_tables/{module}_commF.txt"
    comT=open(commTablepath, "w")
    print(f"Community\tTarget(unknown)\tAnnotated\tFunctions", file=comT) #header
    
    #Communities
    louvain_comms=nx.community.louvain_communities(nw, weight="correlation", seed=42)
    num_communities=len(louvain_comms)
    print(f"Number of communities: {num_communities}")
    
    #Draw each community
    for i, comm in enumerate(louvain_comms):
        
        #ignore community restrictions, obtain functions for all genes
        #functions dictionary for the community
        fd=gan.main(module, comm, 10, 10) 
        #save info to table for easy parsing, also get known and unknown genes
        [uk_G, k_G]=tableWriter(i, fd, comT)
        
        if len(comm)>community_limit:
            #if i==5: #look into a particular community
                subgraph=nw.subgraph(comm)
                s_pos=nx.spring_layout(subgraph, seed=42)
                
                plt.figure(figsize=(width_in, height_in), dpi=dpi)
                plt.title(f"Community {i}: {len(comm)} genes")
                nx.draw(subgraph, pos=s_pos, with_labels=False, node_size=100, edge_color="gray")
                nx.draw_networkx_nodes(subgraph, pos=s_pos, node_size=100, linewidths=1, edgecolors="black")
                
                #color unknown genes black
                unknown_genes={node: node for node in subgraph.nodes if node in uk_G}
                nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=unknown_genes, node_color="black", node_size=100)
                
                #select and color closely related genes
                #nestList_cg=[val[0] for val in fd.values()] #get only corr genes
                #corrGen=sum(nestList_cg, []) #unnest the corr genes list
                #k_genes={node: node for node in subgraph.nodes if node in corrGen} #find out closely correlated genes for a gene(find out if a gene in a community has no close correls)
                #nx.draw_networkx_nodes(subgraph, pos=s_pos, nodelist=k_genes, node_color=f"red", node_size=100) #view closerly correlated genes
                
                known_patch=matplotlib.patches.Patch(label=f'Known genes: {len(k_G)}')
                unknown_patch=matplotlib.patches.Patch(color='black', label=f'Unknown genes: {len(uk_G)}')
                plt.legend(handles=[known_patch, unknown_patch], loc="upper right")

                plt.savefig(f"./genes/graphFigs/{module}_{i}.png", dpi=dpi, bbox_inches="tight") #save plot
                #plt.show()
                plt.close()
                
                #grep C1 turquoise_commF.txt |cut -f4 -d$'\t'|cut -d'.' -f1|sort|uniq -c|sort -n

    comT.close()
    

    #Draw the commmunity separated network
    #Create color mapping
    patchList=[]
    cmap=matplotlib.cm.hsv if num_communities > 10 else matplotlib.cm.tab10
    #assign colors to nodes based on community
    node_colors=[]
    patchList=[]
    for node in nw.nodes():
        for i, community in enumerate(louvain_comms):
            if node in community:
                node_colors.append(cmap(i/max(num_communities-1, 1)))
    
    for i, comm in enumerate(louvain_comms):
        c_patch=matplotlib.patches.Patch(label=f'Community {i}: {len(comm)}', color=cmap(i/max(num_communities-1, 1)))
        patchList.append(c_patch)


    #draw
    plt.figure(figsize=(dimlist[0], dimlist[1]), dpi=dimlist[2])
    nw_pos=nx.spring_layout(nw, seed=42)
    nx.draw_networkx_edges(nw, pos=nw_pos, edge_color="gray", alpha=0.5)
    nx.draw_networkx_nodes(nw, pos=nw_pos, node_color=node_colors, node_size=40, linewidths=1, edgecolors="black")
    
    plt.axis('off')
    plt.title(f"{module.title()} module: {num_communities} communities")
    plt.legend(handles=patchList, loc="upper right")
    plt.savefig(f"./genes/graphFigs/{module}_communities.png", dpi=dpi, bbox_inches="tight") #save plot
    plt.show()

if __name__ == "__main__":
    module="violet" #black, brown
    main(module, 0.6)