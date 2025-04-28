#!/bin/python3

import subprocess

def mtd(c1, c2):
    l1=c1.split() #c1 is sample PE75_flag_4
    l2=c2.split() #c2 is specie PE75
    print("Sample list:",l1)
    print("")
    print("Specie list:",l2)
    
    for i in range(len(l1)):
        if l2[i] not in l1[i]: #c2[i] is specie (short) and c1[i] is sample name, contain specie
            #some species are written with a number in species column, but without number in sample name
            
            word1=l1[i].split("_")[0] #as sample format is DK_26_1, and species is DK150, we take sample DK
            
            if word1 not in l2[i]: #and compare against specie (is DK in DK150?)
                #there are some that have it on the second sample slot B1_EP_1 and species EP1
                word2=l1[i].split("_")[1]
                
                if word2 not in l2[i]: #so we do the same
                
                    print(l1[i], '\t', l2[i]) # sample DK_25_2 	 specie B73
            
    
    return 0

if __name__ == "__main__":
     
    #Metadata is /home/janek/Desktop/Intership/work/data/sample.tissue.correct.cluadj.txt
    c1=subprocess.check_output("cut -f1 /home/janek/Desktop/Intership/work/data/sample.tissue.correct.cluadj.txt", shell=True, text=True)
    c2=subprocess.check_output("cut -f2 /home/janek/Desktop/Intership/work/data/sample.tissue.correct.cluadj.txt", shell=True, text=True)
    
    mtd(c1, c2)