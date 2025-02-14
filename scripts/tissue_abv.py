#!/bin/python3

import sys

def abv_mapping(filepath, endpath):
    tissue_abv={ #dict defuinition, code: abreviation
        "19":"Se10",
        "20":"Se15",
        "21":"Se20",
        "22":"Se25",
        "23":"Se30",
        "24":"Se40M",
        "25":"TaIm",
        "26":"TaMe",
        "27":"CobIm",
        "28":"CobPp",
        "29":"Silk",
        "cob":"Cob",
        "CO":"Coleo",
        "CR":"RoCr",
        "ELO":"Le5El",
        "EZ":"RoPrEl",
        "flag":"LeFl",
        "INT":"FElIn",
        "L1BL":"Le1",
        "L3BL":"Le3Bl",
        "L3SH":"Le3Sh",
        "L5BL":"Le5",
        "L8BL":"Le8",
        "MAT":"Le5M",
        "MER":"Le5Mer",
        "M":"Msphy",
        "MZ":"RoPrMers",
        "PR":"RoPr",
        "RZ":"RoHaZ",
        "SR":"RoSe"
        }
    
    with open(filepath, "r") as file, open(endpath, "w") as endF:
        print("Writing table with decodified abreviations...")
        
        for line in file:
            if "V2" in line: #1rst line has "" at the beggining, marks empty name, replace with - to not lose it
                line=line.replace('""', "-")
            
            listed_line=line.replace('"', "").strip().split("	") #transform lines into lists, separate by tab to keep last column
            #files present some columns with "", remove them for better data clarity
            
            current_item=listed_line[3] #slot 3 in the list (the list is the current line) belongs to tissue abreviation
            
            if current_item in tissue_abv: #just in case
                listed_line[3]=tissue_abv[listed_line[3]] #replace the code for the abreviation
                
            print('\t'.join(listed_line), file=endF) #print each line (columns separated by tab) to the file
    
    print(f"DONE: Codes replaced by abreviations can be found in {endpath}")
    
    return 

def main():
    
    filepath=sys.argv[1]
    endpath=sys.argv[2]
    
    abv_mapping(filepath, endpath)

if __name__ == "__main__":
    main()