#!/usr/bin/env python3

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
        #set Z
        g=g.str.replace("z", "Z")
        #find _ and remove suffix
        g=g.str.split("_").str[0]
        return g
    
    def get_types(self): #has to be split in mercator, prot-scriber, swissprot and original description by &, som elines are anomalous and present extra &
        #3-hydroxyacyl-CoA dehydrogenase *(AIM1/MFP2) 
        #adenylylsulfatase *(FHIT) 
        
        #use & prot, & swiss and & original as separators
        return self.gene[3].split(":")