#!/usr/bin/env python3

from gprofiler import GProfiler
import matplotlib
matplotlib.use('gtk3agg') 
import matplotlib.pyplot as plt

#read list frmo file and put into dic by module
def load_genes(filepath):
    geneList={}
    with open(filepath, 'r') as f:
        for line in f:
            line=line.replace('"', "").split()
            if len(line)>1:
                module=line[1]
                gene=line[0]
                if module in geneList.keys():
                    geneList[module]+=[gene]
                else:
                    geneList[module]=[gene]
    return geneList

#Functional enrichment
def functional_analysis(genes, organism):
    gp=GProfiler(return_dataframe=True)
    
    result=gp.profile(organism=organism, query=genes)
    return result

def save_results(df, output_path):
    df.to_csv(output_path, index=False)
    
def genesToGo(total, annotated, module):
    r=total-annotated
    plotParts=[annotated, r]
    col=["darkgreen", "darkred"]
    label=[f"Annotated genes:\n{annotated}", f"Unknown genes:\n{r}"]
    expand=(0.1, 0)
    
    plt.figure(figsize=(8,6))
    plt.title(f"Annotated gene ratio in ME {module}: {round((annotated/total)*100, 4)}%") #title with % annotated
    plt.pie(plotParts, labels=label, explode=expand, colors=col)
    plt.axis('equal')
    plt.show()

#ME big>small
"turquoise", "blue", "grey", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple"
"greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen", "lightyellow"
"royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white", "skyblue"
"saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen", "darkmagenta", "sienna3"
"yellowgreen", "skyblue3", "plum1"

def main():
    module="turquoise"
    geneFile="../coexpression_code/geneModule.txt"
    outFile=f"MEann/f_ann_{module}.csv"
    organism="zmays"
    
    print(f"Module is: {module}")
    
    geneList=load_genes(geneFile)

    geneMElist=geneList[module]
    print(f"Number of genes is {len(geneMElist)}")
    print(geneMElist[:10])
    
    #use gprofiler to annotate all possible genes
    results=functional_analysis(geneMElist, organism)
    print(f"Number of analyzed genes is {len(results)}")
    
    #save results
    if not results.empty:
        #save_results(results, outFile)
        print(f"Results saved to {outFile}")
    else:
        print("No functional enrichment found.")
    
    #Plot analyzed genes
    genesToGo(len(geneMElist), len(results), module)
    

if __name__ == "__main__":
    #cut -d" " -f2 geneModule.txt|sort|uniq -c
    
    main()
    
'''
4926 "turquoise"
4671 "blue"
3978 "grey"
3637 "brown"
2423 "yellow"
2027 "green"
1361 "red"
1266 "black"
1120 "pink"
1062 "magenta"
990 "purple"
923 "greenyellow"
734 "tan"
580 "salmon"
537 "cyan"
515 "midnightblue"
454 "lightcyan"
409 "grey60"
356 "lightgreen"
350 "lightyellow"
344 "royalblue"
301 "darkred"
282 "darkgreen"
271 "darkturquoise"
271 "darkgrey"
268 "orange"
261 "darkorange"
238 "white"
236 "skyblue"
208 "saddlebrown"
153 "steelblue"
118 "paleturquoise"
83 "violet"
76 "darkolivegreen"
64 "darkmagenta"
44 "sienna3"
37 "yellowgreen"
32 "skyblue3"
28 "plum1"
1 "modules"
'''