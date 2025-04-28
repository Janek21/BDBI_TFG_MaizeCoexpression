#!/usr/bin/env python3

import pandas as pd
import numpy as np

#read ID column, remove _p..., then make dictionary for modules, comparing functions/module and functions/total

class Annotater(object):
    
    def __init__(self, annFile, MEfile):
        
        #creates global ann_data and MEg_dict
        self.data_prep(annFile, MEfile)
        
        print(self.MEg_dic) #blue: genes: g1, g2, g3
        print(self.ann_data)
        
        #self.ann_data.to_csv("./fff.txt", sep='\t')


        self.data_calculation()
        #print(self.MEg_dic)
        
        
    
    def data_prep(self, annFile, MEfile):
        #read file
        ann_data=pd.read_table(annFile, sep="	")
        #remove ' from all the table, for better processing later
        for col in ann_data.columns:
            ann_data[col]=ann_data[col].str.replace("'", "") 
            
        #fix ID column (include removing ')
        ann_data["IDENTIFIER"]=self.id_rectifier(ann_data["IDENTIFIER"])
        #globalize the variable
        self.ann_data=ann_data
        
        #Set up ME-gene link
        #read gene-modules file
        ME_data=pd.read_table(MEfile, sep=" ", names=["genes", "MEs"]) #get 2 columns instead of colnames
        ME_data=ME_data.drop(0) #Remove the first row (old column names)
        
        #create a list and a (global)dictionary to relate modules and genes
        [MElist, self.MEg_dic]=self.MEd_genes(ME_data)
        
        #add column for ME to annotation data
        self.ann_data["MEs"]=MElist
        
        
        



    def data_calculation(self):
        #join by function, functions too extensive (8653), most present minmal differentiation
        #Repetitions of each function in all the dataset
        total_count=self.dic_calc(self.ann_data["NAME"]) #f1:10, f2:3
        print(sum(total_count.values())) ##includes all empty cases and all nonannotation cases, valid?
        
        for module in self.MEg_dic.keys(): #add functions reps of all modules to dictionary
            self.MEd_functions(module)
        #dic is: plum1: {genes:g1,g2}, {functions: {f1:1}, {f2:2}}
        #print(self.MEg_dic["plum1"]["functions"].values())
        print(any("genes" in d for d in self.MEg_dic.values()))
        
        ## How specific do functions have to be (first, 2nd or all orders (.))
        m=self.ratio_calculator("Cell division.cell cycle organisation.metaphase to anaphase transition.transcription factor", self.MEg_dic, "plum1")
        print(m)
        t=self.ratio_calculator("Cell division.cell cycle organisation.metaphase to anaphase transition.transcription factor", total_count)
        print(t)
        print("shorttt")
        t=self.ratio_calculator("Cell division.cell cycle", total_count)
        print(t)
        
        return


    #########################

    def id_rectifier(self, id_col): #set Z and remove suffix to be able to find it in ME comparison
        #set Z
        id_col=id_col.str.replace("z", "Z")
        #find _ and remove suffix
        id_col=id_col.str.split("_").str[0]
        return(id_col)

    def MEd_genes(self, ME_data):
        MElist=[]
        MEg_dic={}
        
        #create the list to add as a column
        for gene in self.ann_data["IDENTIFIER"]: #search all genes in order
            if gene in list(ME_data["genes"]): #if gene exists (with ME) #transform IDs into list so that they can be found
                module = ME_data.loc[ME_data["genes"] == gene, "MEs"].iloc[0]
                
                MElist.append(module)
                
                #use the parsing to create a dictionary
                if module in MEg_dic.keys(): #if the module key is created add the gene
                    MEg_dic[module]["genes"]+=[gene]
                else: #create the new module
                    MEg_dic[module]={"genes":[gene]}
                    
            else: #if ann_ID is a blank entry, add a blank entry as well
                MElist.append("")
                
        return [MElist, MEg_dic]



    def dic_calc(self, func_list):
        #count repetitions of each function
        func_dic={}
        
        for func in func_list: #change to be a single ME or other options
            if func in func_dic.keys():
                func_dic[func]+=1
            else:
                func_dic[func]=1
                
        return func_dic
    
    def MEd_functions(self, module="NULL"): #given a module it introduces all of its function repetitions into MEg_dic
        
        if module=="NULL": #detect errors
            raise ValueError("Missing module to calculate its function repetitions")
        
        #get functions(NAME column) from all genes in ME="blue" (blue in ME>get NAME)
        currentME_functions=self.ann_data.loc[self.ann_data["MEs"] == module, "NAME"]
        #find the function repetitions of the current module
        currentME_count=self.dic_calc(currentME_functions)
        #add it to the dictionary (all function reps) under key Functions
        self.MEg_dic[module]["functions"]=currentME_count
        
    def ratio_calculator(self, c_function, in_dic, module="NULL"): #calculate f1 ratio 
        #detect if total or module dictionary
        if module!="NULL": #module is no chosen, is it an error?
            if any("genes" in d for d in in_dic.values()): #no module+genes as key -> module dictionary
                total=sum(in_dic[module]["functions"].values()) #total functions in mocule
                c_function_reps=in_dic[module]["functions"][c_function] #function reps in current module
                print("ratCalc is ", c_function_reps, total)
                r=c_function_reps/total
                
            else: #no genes as key but no module chosen means an error
                raise ValueError("Missing module to calculate ratio")
        
        else: #the total dictionary is chosen
            r=in_dic[c_function]/sum(in_dic.values())
            print("ratCalc total is ", in_dic[c_function], sum(in_dic.values()))
        
        
        return r
    
    def ratio_comparator(self, f1): 
        #compare f1 ratio from total_counts and self.MEg_dic[module]["Functions"]
        #total=n of genes(or rows), case=cases with this function- total=n of genes in the module, case=reps of this function in the module
        
        #get module ratio
        
        #get totals ratio
        
        #statistical test
        
        return

def main():
    annFile="../data/annotations/b73.mercator.v4.7.txt"
    MEfile="geneModule.txt"
    
    annoChecker=Annotater(annFile, MEfile)
    #data_prep(annFile, MEfile)


if __name__ == "__main__":
    
    main()