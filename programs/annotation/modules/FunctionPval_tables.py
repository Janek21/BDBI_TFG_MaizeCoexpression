#!/usr/bin/env python3

#gets mercator annotation, removes unimportant genes, writes tables of functions+pvalue in each module (a table for each annotation type)
#Step 1

from FileRead_class import Gene
import subprocess
import scipy.stats as sps

#@ filter out cases where no gene is present, there is only a function
#@ filter unassigned/unnanotated genes
#@ remove enzyme classifications
#@ fischer test -> p value > order functions in each module by pvalue, use top 3 significant
#@ gene concretion

def listmaker(genlist, depth, module=None, MEfile=None): #gets reps of all functions (of chosen depth)
    dic_functions={}
    for g in genlist: #iterate through all genes
        
        if "not assigned" in g.get_annotation(): #skip genes with no function
            continue
        
        ME=""
        if module!=None: #if module is not given, the function calculates global, so ME is not needed
            ME=g.get_ME(MEfile) #get ME
        
        if module==None or module==ME: #module==None when no module is specified(calculate for all), module==ME when module is of interest and ME matches
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

def extramaker(genlist, types, module=None, MEfile=None): #gets reps of all functions of the chosen annotation type
    dic_functions={}
    for g in genlist: #iterate through all genes
        
        ann_types=g.get_types()
        
        #get the annotation type of interest
        if ann_types==None: #if its None it wont work for the following "if" results
            continue #None cases are empty cases, where that description is not found, we skip them
        elif types=="mercator":
            ann_types=ann_types[0].strip()
        elif types=="prot-scriber":
            ann_types=ann_types[1].strip()
        elif types=="swissprot":
            ann_types=ann_types[2].strip()
        elif types=="original":
            ann_types=ann_types[3].strip()
        else:
            raise ValueError("The possible annotation types are: mercator, prot-scriber, swissprot and original")
        
        if "no annotation" in ann_types or "not classified" in ann_types or "none" in ann_types: #"no annotation" > swissprot, prot-scriber | "not classified" > mercator | "none" > original description
            continue #skip genes with no function
        
        ME=""
        if module!=None: #if module is not given, the function calculates global, so ME is not needed
            ME=g.get_ME(MEfile) #get ME #as this process is slow, only do it when it is needed(this is why is the last part of the function)
        
        if module==None or module==ME: #if ME coresponds with target module, or the calculation is for global
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

def fisher_calc(total_dic, module_dic, ann_func):
    alpha=0.05 
    funcME=module_dic[ann_func] #"blue"fA
    nonfuncME=sum(module_dic.values())-funcME #"blue"nonfA=all blue-"blue"fA
    
    funcRest=total_dic[ann_func]-funcME #"nonblue"fA=all fA-"blue"fA
    nonfuncRest=sum(total_dic.values())-sum(module_dic.values())-funcRest #"nonblue"nonfA=all-all blue-"nonblue"fA
    
    fishertable=[[funcME, funcRest], [nonfuncME, nonfuncRest]]
    pval=sps.fisher_exact(fishertable).pvalue
    pval=float(pval)
    
    #H0=The annotated function is not relevant in this module
    #H1=The annotated function is relevant in the module
    
    if alpha>pval:#if pvalue is smaller than alpha, we reject the null (True that the function is relevant)
        return [True, pval]
    
    return [False, pval]


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
    
    MElist=moduleset(modulefile)
    typePossiblities=["mercator", "prot-scriber", "swissprot"] #original is always None(cut -f4 -d$'\t' b73.mercator.v4.7.txt|rev|cut -c 1-27|rev|sort|uniq), so we ignore it for now
    
    
    print(MElist)
    
    for i in range(3):#loop trhough all depths and types
        depth=i+1
        tp=typePossiblities[i]
        print(tp, depth)
        
        fileD=open(f"./annotation/modules/Pres/d_{depth}_sheet.txt", "a") #a file for each depth type
        fileT=open(f"./annotation/modules/Pres/anT_{tp}_sheet.txt", "a") #a file for each annotation type
        
        print(f"Function\tTotalCounts\tModuleCounts\tPvalue\tRelevancy\tModule", file=fileD)
        print(f"Function\tTotalCounts\tModuleCounts\tPvalue\tRelevancy\tModule", file=fileT)
        
        annall_dic=listmaker(genlist, depth) #dict of all function:reps of general annotation
        typesall_dic=extramaker(genlist, tp) #dict of all function:reps from a concrete annotation type
        for m in MElist:
            print(m)
            
            annME_dic=listmaker(genlist, depth, m, modulefile) #get dict of all functions:reps from a module
            for ann_function in annME_dic.keys(): #check which out of all annotated functions are relevant
                res=fisher_calc(annall_dic, annME_dic, ann_function)
                print(f"{ann_function}\t{annall_dic[ann_function]}\t{annME_dic[ann_function]}\t{res[1]}\t{res[0]}\t{m}", file=fileD)
            
            typeME_dic=extramaker(genlist, tp, m, modulefile) #get dict of all functions:reps from a module (from the concrete type of annotation))
            for type_function in typeME_dic.keys(): #check which out of all type annotated functions from the chosen annotation type are relevant
                res=fisher_calc(typesall_dic, typeME_dic, type_function)
                type_function_name=type_function.split(": ")[1] #remove prot-scriber: and others
                print(f"{type_function_name}\t{typesall_dic[type_function]}\t{typeME_dic[type_function]}\t{res[1]}\t{res[0]}\t{m}", file=fileT)
                
        fileD.close()
        fileT.close()

    
if __name__ == "__main__":
    main()