#!/usr/bin/env python3

from gene_class import annotated_function, Gene
import matplotlib.pyplot as plt

def notann_plotter(genlist, enzyme, nann):
    
    nass=0
    for g in genlist:
        if "not assigned" in g.get_annotation(): #skip genes with no function
            nass+=1
    

    entryTypes=[len(genlist)-nass, enzyme, nann, nass]
    total=sum(entryTypes)
    ratios=[round((entryTypes[0]/total)*100, 4), round((entryTypes[1]/total)*100, 4), round((entryTypes[2]/total)*100, 4), round((entryTypes[3]/total)*100, 4)]
    entryLabs=[f"Annotated: {ratios[0]}%", f"Enzymes: {ratios[1]}%", f"Not annotated: {ratios[2]}%", f"Not assigned: {ratios[3]}%"]
    print(entryTypes)
    expand=(0, 0.1, 0, 0)
    plt.pie(entryTypes, labels=entryLabs, explode=expand)
    plt.show()
    
    #return

def main():
    
    genlist=[]
    enzyme=0
    nann=0
    with open("../data/annotation/b73.mercator.v4.7.txt") as f: #open file
        for line in f:
            g=Gene(line)
            
            if "NAME" in g.gene: #filter out header, first check, if not following give errors
                continue
            if g.get_id()==None: #Filter entries with no gene (function headers for tree building), early, as it can give subsequent errors
                continue
            if any("enzyme" in x.lower() for x in g.get_annotation()): #enzymes dont count for functions, but have to be counted
                enzyme+=1
                continue
            if "not annotated" in g.get_annotation(): #dont include completely unnanotated genes, but count them
                nann+=1
                continue
            
            genlist.append(g) #if none of the earlier apply, add to list of gene objects
        
        notann_plotter(genlist, enzyme, nann)
            
    typePossiblities=["mercator", "prot-scriber", "swissprot"] #original is always None(cut -f4 -d$'\t' b73.mercator.v4.7.txt|rev|cut -c 1-27|rev|sort|uniq), so we ignore it for now


if __name__ == "__main__":
    main()