#!/usr/bin/env python3

#Given a funcion and a ME outputs all the corresponding genes
#Step 4
#@get genes only from function??

from modules.FileRead_class import Gene
import pandas as pd

def MEfu(genlist, ME, annFu, MEfile, depth=1): #list of genes with ME and Function
    outGlist=[]
    for g in genlist:
        if g.get_annotation(depth)[0]==annFu:
            if g.get_ME(MEfile)==ME:
                outGlist.append(g.get_id())
    
    return outGlist


def main():
    modulefile="geneModule.txt"
    genlist=[]
    with open("../data/annotation/b73.mercator.v4.7.txt") as f: #open file
        for line in f:
            g=Gene(line)
            
            if "NAME" in g.gene: #filter out header, first check, if not following give errors
                continue
            if g.get_id()==None: #Filter entries with no gene (function headers for tree building), early, as it can give subsequent errors
                continue
            if any("enzyme" in x.lower() for x in g.get_annotation()): ##would a simple "Enzyme classification" in g.get_annotation() be enough? or eliminate enzyme in further depth as well #or go fo all in g.gene
                continue
            if "not annotated" in g.get_annotation(): #dont include completely unnanotated genes (present empty types as well)
                continue
            
            genlist.append(g) #if none of the earlier apply, add to list of gene objects
        
    resGen=MEfu(genlist, "turquoise", "Photosynthesis", modulefile, 1)
    
    with open("../data/annotation/genes_all.txt") as file:
        totGen=file.readlines()
    
    ff=open("./ff.txt", "a")
    matchGen=[]
    for line in totGen:
        line=line.strip().split("\t")
        if line[0] in resGen:
            relevant=[line[0]]+line[3:5]
            matchGen.append(relevant)
            print("\t".join(relevant), file=ff)
    ff.close()
    print(matchGen)
    
    
if __name__ == "__main__":
    main()