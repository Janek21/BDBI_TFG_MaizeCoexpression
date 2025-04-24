#!/usr/bin/env python3

import pandas as pd

#read ID column, remove _p..., then make dictionary for modules, comparing functions/module and functions/total

#class Annotater(object):?
def data_prep(annFile, MEfile):
    #read file
    ann_data=pd.read_table(annFile, sep="	")
    #fix ID column (include removing ')
    ann_data["IDENTIFIER"]=id_rectifier(ann_data["IDENTIFIER"].str.replace("'", ""))
    #remove ' from BINCODE, fpr better processing if needed
    ann_data["BINCODE"]=ann_data["BINCODE"].str.replace("'", "")
    
    #Set up ME-gene link
    #read gene-modules file
    ME_data=pd.read_table(MEfile, sep=" ", names=["genes", "MEs"]) #get 2 columns instead of colnames
    ME_data=ME_data.drop(0) #Remove the first row (old column names)
    
    #create a list and a dictionary to relate modules and genes
    [MElist, MEg_dic]=ME_joiner(ann_data, ME_data)
    
    #add column for ME to annotation data
    ann_data["MEs"]=MElist
    print(ann_data)



def data_calculation():
    #join by function, functions too extensive (8653), most present minmal differentiation
    func_dic=dic_calc(ann_data)
    
    return


#########################

def id_rectifier(id_col): #set Z and remove suffix to be able to find it in ME comparison
    #set Z
    id_col=id_col.str.replace("z", "Z")
    #find _ and remove suffix (readd ' beacuse it is erased as well)
    id_col=id_col.str.split("_").str[0]
    return(id_col)

def ME_joiner(ann_data, ME_data):
    MElist=[]
    MEg_dic={}
    
    #create the list to add as a column
    for gene in ann_data["IDENTIFIER"]: #search all genes in order
        if gene in list(ME_data["genes"]): #if gene exists (with ME) #transform IDs into list so that they can be found
            module = ME_data.loc[ME_data["genes"] == gene, "MEs"].iloc[0]
            
            MElist.append(module)
            
            #use the parsing to create a dictionary
            if module in MEg_dic.keys(): #if the module key is created add the gene
                MEg_dic[module]+=[gene]
            else: #create the new module
                MEg_dic[module]=[gene]
                
        else: #if ann_ID is a blank entry, add a blank entry as well
            MElist.append("")
            
    return [MElist, MEg_dic]



def dic_calc(ann_data):
    #count numbers of each function
    func_dic={}
    for func in ann_data["NAME"]:
        
        if func in func_dic.keys():
            func_dic[func]+=1
        else:
            func_dic[func]=1
    return func_dic

def main():
    annFile="../data/annotations/b73.mercator.v4.7.txt"
    MEfile="geneModule.txt"
    data_prep(annFile, MEfile)


if __name__ == "__main__":
    
    main()
    
    
    
#ann_ID=ann_data["IDENTIFIER"].str.replace("'", "") #remove all ', genes cant be found if ' present
#ME_joiner(ann_ID, ME_data, 0, [], {})
def ME_joiner_recursive(ann_ID, ME_data, l, MElist, MEg_dic):
    
    #create the list to add as a column
    if l<len(ann_ID):
        gene=ann_ID[l]
        print(l)
        
        if gene in list(ME_data["genes"]): #if gene exists (with ME) #transform IDs into list so that they can be found
            module = ME_data.loc[ME_data["genes"] == gene, "MEs"].iloc[0]
            
            #use the parsing to create a dictionary
            if module in MEg_dic.keys(): #if the module key is created add the gene
                MEg_dic[module]+=[gene]
            else: #create the new module
                MEg_dic[module]=[gene]
            
            return ME_joiner(ann_ID, ME_data, l+1, MElist+[module], MEg_dic)
            
                
        else: #if ann_ID is a blank entry, add a blank entry as well
            return ME_joiner(ann_ID, ME_data, l+1, MElist+[""], MEg_dic)
    
    return MElist, MEg_dic