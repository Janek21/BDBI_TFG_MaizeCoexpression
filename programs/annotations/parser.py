#!/usr/bin/env python3

from gene_class import Gene

with open("../data/annotations/b73.mercator.v4.7.txt") as f:
    tt=[]
    for line in f:
        tt.append(line)

print(Gene(tt[5]).get_types())