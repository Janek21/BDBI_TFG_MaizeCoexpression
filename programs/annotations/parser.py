#!/usr/bin/env python3

from gene_class import Gene
import random, subprocess

#@ filter unassigned/unnanotated genes [yes/no]
#@ filter out cases where no gene is present, there is only a function

def listmaker(genlist, depth, module=None, MEfile=None):
    dic_functions={}
    for g in genlist: #iterate through all genes
        
        ME=""
        if module!=None: #if module is not given, the function calculates global, so ME is not needed
            ME=g.get_ME(MEfile) #get ME (brown or None)
        
        if module==None or module==ME: #module==None when no module is specifies(calculate for all), module==ME when module is of interest and ME matches
            functs=g.get_annotation(depth) #get [depth] levels of annotation
            
            if depth>1 and depth==len(functs): #if depth > 1 it will be joint (2+ segment annotations are treated as 1 str element)
                functs=".".join(functs)
            else:  #if its 1 element, get it turned to string
                functs=functs[0]
            
            if functs in dic_functions.keys(): #if the function has already been found add 1 to counter, if its new, set it to 1
                dic_functions[functs]+=1
            else:
                dic_functions[functs]=1
        
    return dic_functions #return reps of all functions (all or all from one ME)

def extramaker(genlist, types, module=None, MEfile=None):
    dic_functions={}
    for g in genlist: #iterate through all genes
        
        ME=""
        if module!=None: #if module is not given, the function calculates global, so ME is not needed
            ME=g.get_ME(MEfile) #get ME (brown or None) 
        
        if module==None or module==ME:
            ann_types=g.get_types()
            
            #get the annotation type of interest
            if ann_types==None: #if its None it wont work for the following if results
                ann_types="None" #None cases are empty cases, where that description is not found
            elif types=="mercator":
                ann_types=ann_types[0].strip()
            elif types=="prot-scriber":
                ann_types=ann_types[1].strip()
            elif types=="swissprot":
                ann_types=ann_types[2].strip()
            elif types=="original":
                ann_types=ann_types[3].strip()
            else:
                return "The possible annotation types are: mercator, prot-scriber, swissprot and original"
            
            #Count the repetitions of all annotation types
            if ann_types in dic_functions.keys():
                dic_functions[ann_types]+=1
            else:
                dic_functions[ann_types]=1
            
    return dic_functions

def moduleset(MEfile):
    modulecolumn=subprocess.run(["cut", "-d", " ", "-f2", MEfile], capture_output=True) #get the 2nd column of the file, where the mocules are listed
    MEset=modulecolumn.stdout.decode().replace('"', "").split() #deode the column, remove all the " and turn it into a list
    MEset=set(MEset[1:]) #transform into a set to remove repeated elements and remove the first one as its the column title
    return MEset

def ratio(d_func, function):
    
    #d_func=nann_rm(d_func) ##remove unannotated?
    
    f_case=d_func[function] #repetitions of the specified function
    total=sum(d_func.values()) #total amount to calculate ratio 
    r=f_case/total
    return r

##is it needed? should we keep or not non assgned/non annotated genes
def nann_rm(d_func): #remove not assigned(not annotated and annotated) elements
    for key in list(d_func.keys()):
        if "not assigned" in key:
            d_func.pop(key)
    return d_func

def main():
    
    modulefile="geneModule.txt"
    genlist=[]
    with open("../data/annotations/b73.mercator.v4.7.txt") as f: #open file
        for line in f:
            g=Gene(line)
            if not g.get_id()==None: ##Filter None ID cases? (cases that only present function) (2000+ nonexsitent genes are out)
                genlist.append(g) #make a list of gene objects
    
    if "NAME" in genlist[0].gene: #if the first line is the header, remove it
        genlist.pop(0)

    ## For loop each block to analyze all functions in 1 module (blocks are normal functions and extra functions)
    ## Add a top for loop to iterate trhough all modules as well and write it all to file
    ## Another extra loop can be added for the depth and annType as well (a different file for each)
    
    ## Module of each function is written in the row (or maybe a directory+infilename for each module)
    ## depth/annType level is written in file name
    
    typePossiblities=["mercator", "prot-scriber", "swissprot", "original"]
    
    for i in range(4):
        anntype=typePossiblities[i]
        depth=i #ignore when 0
    #create dictionaries for the repetitions of functions and annotations belonging to a specific type
    d_ann_tt=listmaker(genlist, 3)
    d_annT_tt=extramaker(genlist, "prot-scriber")
    
    MEset=moduleset(modulefile) #get a list of all the possible modules, no repetitions
    for m in MEset:
        
        #create dictionaries for the repetitions of functions and annotation types for a specific module
        d_ann_me=listmaker(genlist, 3, module=m, MEfile=modulefile)
        d_annT_me=extramaker(genlist, "prot-scriber", module=m, MEfile=modulefile)
        
        
        #Calculate ratio, representation and write to file for functions
        with open(f"./annotations/{m}test.txt", "a") as file:
            print(f"Function\tMEratio\tMEtotal\tTratio\tTtotal\tME_representation", file=file) #set header
            for function_i in list(d_ann_me.keys()): #calculate the representation of all functions in the module
                print("Function is:", function_i)
                
                r_me=ratio(d_ann_me, function_i) #calculate the ratio of teh function in the module
                r_tt=ratio(d_ann_tt, function_i) #calculate the ratio of the function in total
                
                #write file
                print(f"{function_i}\t{r_me}\t{d_ann_me[function_i]}\t{r_tt}\t{d_ann_tt[function_i]}\t{r_me>=r_tt}", file=file)
        
        
        #Calculate ratio, representation and write to file for different annotation types
        with open(f"./annotations/T{m}test.txt", "a") as fileT:
            print(f"Function\tMEratio\tMEtotal\tTratio\tTtotal\tME_representation", file=fileT) #set header
            for functionT_i in list(d_annT_me.keys()): #calculate the representation of all the extra annotations of 1 type in the module
                print("Function is:", functionT_i)
                
                rT_me=ratio(d_annT_me, functionT_i)
                rT_tt=ratio(d_annT_tt, functionT_i)
                
                #Write file
                functionName=functionT_i.split(": ")[1] #Dont get the annotation source in the table
                print(f"{functionName}\t{rT_me}\t{d_annT_me[functionT_i]}\t{rT_tt}\t{d_annT_tt[functionT_i]}\t{rT_me>=rT_tt}", file=fileT)
    
    
    '''sorted_dict={}
    for key in sorted(d_annT_me, key=d_annT_me.get):
        sorted_dict[key]=d_annT_me[key]
    print(sorted_dict)'''
    
if __name__ == "__main__":
    main()