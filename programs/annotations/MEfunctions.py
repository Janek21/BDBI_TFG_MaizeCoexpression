#!/usr/bin/env python3

from gene_class import Gene, annotated_function

def ME_group(annlist):
    overrepME_dic={}
    for af in annlist:
        relevance=af.get_rel() #Relevant in module True/False
        if relevance: #if relevant, they are valid as module functions
            ME=af.get_ME() #get ME and function
            funct=af.get_function()
            pval=af.get_pval()
            
            if ME in overrepME_dic.keys(): #set up a dict for the representative functions of all modules
                overrepME_dic[ME].append([funct, pval])
            else:
                overrepME_dic[ME]=[[funct, pval]]
            
    return overrepME_dic

def diKeyer(representationD, ME):
    ls=representationD[ME]
    ls.sort(key=lambda x: x[1])
    print(ME)
    print(*ls[:5])
    print("")
    

def main():
    filepath="./annotations/Pres/anT_mercator_sheet.txt"
    filepath="./annotations/Pres/d_1_sheet.txt"
    annlist=[]
    with open(filepath) as fp:
        for line in fp:
            l=annotated_function(line)
            annlist.append(l)
    
    if "TotalCounts" in annlist[0].entry:
        annlist.pop(0)
    
    repME=ME_group(annlist) #"brown"
    for ME in repME.keys():
        diKeyer(repME, ME)


if __name__ == "__main__":
    main()