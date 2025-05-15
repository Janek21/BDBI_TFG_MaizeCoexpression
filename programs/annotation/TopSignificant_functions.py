#!/usr/bin/env python3

#writes a table of most significant functions+module in all annotation types

from gene_class import annotated_function

def ME_group(annlist):
    overrepME_dic={}
    for af in annlist:
        relevance=af.get_rel() #Relevant in module True/False
        if relevance: #if relevant, they are valid as module functions
            ME=af.get_ME() #get ME, function and pval
            funct=af.get_function()
            pval=af.get_pval()
            
            if ME in overrepME_dic.keys(): #set up a dict for the representative functions of all modules
                overrepME_dic[ME].append([funct, pval]) #dict contains ME:function, pval
            else:
                overrepME_dic[ME]=[[funct, pval]]
            
    return overrepME_dic

def diKeyer(representationD, ME, toplevel=3): #orders dictionary values by significance(pval), and geths the top x most significant
    ls=representationD[ME] #get [functions:pval] for chosen module
    ls.sort(key=lambda x: x[1]) #order
    return(ls[:toplevel])
    

def main():
    filepath="./annotation/Pres/anT_swissprot_sheet.txt"
    #filepath="./annotation/Pres/d_3_sheet.txt"
    functiontype="swissprot" #d1,d2,d3, mercator, prot-scriber, swissprot
    annlist=[]
    siglevel=3 #how many levels of significant functions (top 3, etc)
    with open(filepath) as fp:
        for line in fp:
            l=annotated_function(line)
            annlist.append(l)
    
    if "TotalCounts" in annlist[0].entry: #remove header
        annlist.pop(0)
    
    with open("./annotation/Pres/database.txt", "a") as datafile: #create file to write results
        readfile=open("./annotation/Pres/database.txt", "r") #read file to check if there is header
        fullfile=readfile.read(1)
        if not fullfile: #if the file is empty(first time writing into it)
            sigFunc=[f"SFunction_{x+1}" for x in range(siglevel)]#create as many function names as significant functions
            functionsHeader='\t'.join(sigFunc)
            print(f"Function_origin\tModule\t{functionsHeader}", file=datafile) #write header
        readfile.close()
        
        repME=ME_group(annlist) #dic of ME:[[functions:pval],etc]
        for ME in repME.keys():
            topF_list=diKeyer(repME, ME, toplevel=siglevel)#takes list of [func, pvalue] for an ME and output top 3 most significant
            nameTopF=[x[0] for x in topF_list] #get names of most significant functions
            r='\t'.join(nameTopF) #space functions to be able to identify in table
            print(f"{functiontype}\t{ME}\t{r}", file=datafile)


if __name__ == "__main__":
    main()