#!/usr/bin/env python3

from gene_class import Gene, annotated_function

def ME_group(annlist):
    overrepME_dic={}
    for af in annlist:
        rep=af.get_rep() #Overrepresented in module True/False
        if rep: #if overrepresented, they are valid as module functions
            ME=af.get_ME() #get ME and function
            funct=af.get_function()
            r=af.get_MEratio()
            
            if ME in overrepME_dic.keys(): #set up a dict for the representative functions of all modules
                overrepME_dic[ME].append([funct, r])
            else:
                overrepME_dic[ME]=[[funct, r]]
            
    return overrepME_dic

def diKeyer(representationD, ME):
    ls=representationD[ME]
    ls.sort(key=lambda x: x[1], reverse=True)
    print(ME)
    print(*ls[:11])
    print("")
    

def main():
    filepath="./annotations/anT_prot-scriber_sheet.txt"
    filepath="./annotations/d_1_sheet.txt"
    annlist=[]
    with open(filepath) as fp:
        for line in fp:
            l=annotated_function(line)
            if "MEratio" not in l.entry:
                annlist.append(l)
    
    repME=ME_group(annlist) #"brown"
    for ME in repME.keys():
        diKeyer(repME, ME)


if __name__ == "__main__":
    main()