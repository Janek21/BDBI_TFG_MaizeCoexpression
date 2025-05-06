#!/usr/bin/env python3

import subprocess

class Gene(object):
    
    def __init__(self, line):
        line=line.replace("'", "") #remove ' elements
        line=line.strip("\n").split('	') #turn the line into a list for better parsability
        self.gene=line #define a global variable
        
    def get_code(self): #binary code is the first element of line
        return self.gene[0]
    
    def get_annotation(self, depth=None): #obtain the annotation of that entry
        g=self.gene[1]
        g=g.split(".") #different annotation levels, each one a position in a list general>concrete
        
        if depth is None: #return the whole annotation unless otherwise specified
            depth=len(g)
        return g[:depth]
    
    def get_id(self): #get id of gene, and clean it up
        g=self.gene[2]
        
        if g=="": #No id, then return None
            return None
        #set Z
        g=g.replace("z", "Z")
        #find _ and remove suffix
        g=g.split("_")[0]
        return g

    def get_types(self): #has to be split in mercator, prot-scriber, swissprot and original description by &, some lines are anomalous and present extra &
        g=self.gene[3]
        
        if g=="": #if there is no decription of the types annotations, return None
            return None
        #use & prot, & swiss and & original as separators
        mercator_split=g.split(" & prot-scriber: ") # 2 blocks, first is mercator, 2nd is the rest
        protscriber_split=mercator_split[1].split(" & swissprot: ") #from the previous rest block get 2 blocks, 1st one protscriber, 2nd is rest
        swissprot_split=protscriber_split[1].split(" & original description: ") #from the rest block split by original desc., get 2 blocks 1st is swissprot, 2nd is original decription(rest)
        
        #from each of the lists get the first element(concrete annotation type), except original, that is the second element
        #readd the type of annotation
        mercator=mercator_split[0]
        protscriber="prot-scriber: "+protscriber_split[0]
        swissprot="swissprot: "+swissprot_split[0]
        original="original description: "+swissprot_split[1]
        
        return [mercator, protscriber, swissprot, original]
    
    def get_ME(self, MEfile): ##slowest part, for loop is faster than grep?
        g=self.get_id() #get gene ID
        
        if g==None:
            return None
        
        grep_out=subprocess.run(["grep", g, MEfile], capture_output=True) #use grep to find the ID in the gene-module file
        
        if grep_out.returncode==1: #returncode=1 means that the gene was not found in the file
            return None
        
        ME=grep_out.stdout.decode().split() #convert the gene-module line into a list
        ME=ME[1].replace('"', "") #Extract the module
        return ME

  
class annotated_function(object):
    def __init__(self, line):
        line=line.strip("\n").split('\t') #turn the line into a list for better parsability
        self.entry=line #define a global variable
        
    def get_function(self):
        an_f=self.entry[0]
        return an_f
    def get_rep(self):
        an_f=bool(self.entry[5])
        return an_f
    def get_ME(self):
        an_f=self.entry[6]
        return an_f
    
    def get_MEratio(self):
        an_f=self.entry[1]
        return float(an_f)
    def get_MEtotal(self):
        an_f=self.entry[2]
        return int(an_f)
    def get_Tratio(self):
        an_f=self.entry[3]
        return float(an_f)
    def get_Ttotal(self):
        an_f=self.entry[4]
        return int(an_f)
    