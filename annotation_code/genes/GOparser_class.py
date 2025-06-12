#!/usr/bin/env python3

from goatools.obo_parser import GODag
import subprocess

class GOfile(object):
    def __init__(self, line):
        line=line.strip().split("\t")
        self.entry=line
    
    def get_ID(self):
        id=self.entry[1]
        return id
    
    def get_GO(self):
        go_term=self.entry[4]
        return go_term
    
    def get_V5ID(self, translation): #the GAF file is in gene ID V4, our data in V5, use genes_all.txt created dic to translate
        v4ID=self.get_ID()
        if v4ID in translation.keys():
            v5ID=translation[v4ID]
            return v5ID
        else:
            return None

def v4tov5_id(filepath): #dictionary that produces equivalance of v4 to 5 notation and inverse using genes_all.txt
    translation_v4to5={}
    with open(filepath, "r") as file:
        for line in file:
            if "Source" not in line: #no header
                line=line.strip("\n").split("\t")
                v5=line[0]
                v4=line[1]
                translation_v4to5[v4]=v5 #it shouls always be a 1 to 1 correaltion, no 2 old ids for 1 new
    
    return translation_v4to5