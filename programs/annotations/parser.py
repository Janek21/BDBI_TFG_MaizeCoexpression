#!/usr/bin/env python3

from gene_class import Gene

#3887

def listmaker(genlist, depth, module=None, MEfile=None):
    dic_functions={}
    for g in genlist: #iterate through all genes
        
        ME=""
        if module!=None: #if module is not given, the function calculates global, so ME is not needed
            ME=g.get_ME(MEfile) #get ME (brown or None) ##slowest part
        
        if module==None or module==ME: #module==None when no module is specifies(calculate for all), module==ME when module is of interest and ME matches
            functs=g.get_annotation(depth) #get [depth] levels of annotation
            
            if depth>1 and depth==len(functs): #if depth > 1 it will be joint (2+ segment annotations are treated as 1 str element)
                functs=".".join(functs)
            else:  #if its 1 element, get it turned to string
                functs=functs[0]
            ''' ##option to filter out any mismatched depth/true length cases
            elif depth==1:
                functs=functs[0]
            else:
                continue
            '''
            
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
            
            #Count the repetitions of all annottation types
            if ann_types in dic_functions.keys():
                dic_functions[ann_types]+=1
            else:
                dic_functions[ann_types]=1
            
    return dic_functions

def ratio(d_func, function):
    
    d_func=nann_rm(d_func)
    
    f_case=d_func[function] #repetitions of the specified function
    total=sum(d_func.values()) #total amount to calculate ratio ##remove unannotated?
    r=f_case/total
    return r

##is it needed? should we keep or not non assgned/non annotated genes
def nann_rm(d_func): #remove not assigned(not annotated and annotated) elements
    for key in list(d_func.keys()):
        if "not assigned" in key:
            d_func.pop(key)
    return d_func

def main():
    
    genlist=[]
    with open("../data/annotations/b73.mercator.v4.7.txt") as f: #open file
        for line in f:
            if "NAME" not in line: #filter first line (headers)
                genlist.append(Gene(line)) #make a list of gene objects

    d_ann=listmaker(genlist, 3, module="brown", MEfile="geneModule.txt") #create a dictionary of the repetitions of the functions
    
    sorted_dict={}
    for key in sorted(d_ann, key=d_ann.get, reverse=True): #sort by value (big>small)
        sorted_dict[key]=d_ann[key]
    print(sorted_dict)

    #r=ratio(d_ann, "Enzyme classification.EC_2 transferases") #calculate the ratio
    #print(r)
    
    d_annT=extramaker(genlist, "prot-scriber")
    sorted_dict={}
    for key in sorted(d_annT, key=d_annT.get):
        sorted_dict[key]=d_annT[key]
    print(sorted_dict)
    
    g=genlist[23]
    e=g.get_types()
    
if __name__ == "__main__":
    main()