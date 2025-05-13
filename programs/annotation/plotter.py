#!/usr/bin/env python3

from gene_class import annotated_function, Gene
import matplotlib.pyplot as plt
import subprocess
import pandas as pd

## ME and function in axis, do scatterplot, map colors/sizes as well, one of these 2 as pval

def gene_pvaluer(genlist, depth, pvfile, MEfile):
    pvalue_data=pd.read_table(pvfile)
    pvlist=[]
    flist=[]
    for g in genlist:
        ann=g.get_annotation(depth)
        ME=g.get_ME(MEfile)
        ann=".".join(ann)
        rowsAnn=pvalue_data[pvalue_data['Function']==ann]
        matching=rowsAnn[rowsAnn['Module']==ME] #get pvalue of the function of this gene (in the correspondent module)
        relevantmatch=matching[matching["Relevancy"]==True]
        if relevantmatch.empty: #avoid error
            continue
        #print(relevantmatch)
        pval=matching["Pvalue"].values[0]
        pvlist.append(float(pval))
        flist.append(ann)
    return pvlist, flist
        
def scatterplotter(fset, flist, pvlist):
    values = ['a', 'c', 'b', 'd', 'a', 'b', 'c', 'd']
    categories = ['a', 'b', 'c', 'd']
    
    values=flist
    categories=set(flist)
    category_to_y = {cat: i for i, cat in enumerate(categories)}
    m=max(pvlist)
    highlist=[(m-x)*2000 for x in pvlist]
    print(highlist)
    
    y_values = [category_to_y[val] for val in values]
    x_values = list(range(len(values)))


    catplot=plt.scatter(x_values, y_values, s=highlist, c=pvlist, cmap="viridis", alpha=0.2)
    plt.yticks(range(len(categories)), categories)
    plt.grid(True)
    plt.colorbar(catplot, label='Value')
    plt.show()
        

def main():
    dbF="./annotation/Pres/database.txt"
    modulefile="geneModule.txt"
    depth=1
    tann=f"d{depth}"
    pvalues_file=f"./annotation/Pres/{tann.replace('d', 'd_')}_sheet.txt"
    genlist=[]
    with open("../data/annotation/b73.mercator.v4.7.txt") as f:
        #create a list of gene objects (avoid header, empty lines, unnanotated genes and genes involved with Enzymes, as they are not useful for function annotation)
        for line in f:
            g=Gene(line)
            
            if "NAME" in g.gene:
                continue
            if g.get_id()==None:
                continue
            if any("enzyme" in x.lower() for x in g.get_annotation()):
                continue
            if "not annotated" in g.get_annotation():
                continue
            
            genlist.append(g)
        
    checkOut=f"grep {tann} {dbF}|cut -d$'\t' -f3,4"
    output = subprocess.check_output(["bash", "-c", checkOut])
    outrefined=output.decode().replace("\t", "\n").split("\n")
    f_set=set(outrefined)
    print(len(outrefined), len(f_set))
    print(f_set)
    pv_list, f_list=gene_pvaluer(genlist, depth, pvalues_file, modulefile)
    print(pv_list)
    print("AAAAAAAAA")
    scatterplotter(f_set, f_list, pv_list)
    
    ##iterate through each gene (in module X? (then repeat for all modules and plot dide by side)) and test pvalue?
            ##indetad of plotting pvalue per gene, get gene function, and get pvalue in that module
    ##find a way to do it for each tissue as well?

if __name__ == "__main__":
    main()