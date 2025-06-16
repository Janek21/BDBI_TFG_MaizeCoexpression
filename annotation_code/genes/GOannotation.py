#!/usr/bin/env python3

#Given a module(name) and a community(gene list), produce unknown genes in the community, their n closely related genes and the GO terms and functions of these genes

import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

import GOparser_class as gpc

def load_genes(filepath): #create dictionary of genes per module
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

def geneV5_GOterm(genelist, golist, translation): #associate GO terms to V5genes in the module
    geneto_goterm={}
    for go in golist: #iterate through go list IDs as v4 IDs they can be translated into v5 to match in genelist, but v5 IDs to v4 presents missing translations
        id=go.get_V5ID(translation) #IDs in genelist are v5, so use a v5 ID for correspondence
        
        if id in genelist: #check if the gene is of interest
            goTerm=go.get_GO() #get GO terms
            
            if id in geneto_goterm.keys(): #if gene in dictionary, add GOterm
                geneto_goterm[id].add(goTerm)
            else:
                geneto_goterm[id]={goTerm} #GO terms in sets to avoid repetition as gene ID=collection of GO terms
                
    return geneto_goterm #dictionary of gene:[GO terms]
    
def term_function(gene2GO, obo_GOlist): #extract the functions of the GO terms
    gene_function={}
    for gene in gene2GO.keys(): #iterate through the genes
        for go in gene2GO[gene]: #for each GO term belonging to a gene
            
            if go in obo_GOlist and gene in gene_function.keys(): #known term and gene in the dict already
                gene_function[gene]+=[obo_GOlist[go].name]
                
            elif go not in obo_GOlist and gene in gene_function.keys(): #unknown term and gene in the dict already
                gene_function[gene]+=["Unknown"]
                
            elif go in obo_GOlist and gene not in gene_function.keys(): #known term and gene not in the dict already
                gene_function[gene]=[obo_GOlist[go].name]
            
            else: #unknown term+new gene
                gene_function[gene]=["Unknown"]
    
    return gene_function #dictionary of gene:[GO functions]

def unknownGene(genelist, functionslist): #detect genes that have no v4 equivalence, so no go GO terms either, they are unknown
    ukwG=[] #list of unknown genes
    for g in genelist:
        if g not in functionslist: #all genes will be in the dictionary with go functions except unknown ones
            ukwG.append(g)
    return ukwG
    
def closeToGene(unknown_gene, correlation_table, n=5): #gievn an unknown gene ad a correaltion table, detect the n closely correlated genes to the target
    top_cor=correlation_table[unknown_gene].drop(unknown_gene).sort_values(ascending=False).head(n) #get top n correlated genes
    top_corGenes=top_cor.index.tolist() #get the gene names
    return top_corGenes #return the list of closely related genes

def main(module, genelist, n_corrGenes=10, n_sigGO=10): #main function,executes all others, produces the result
    #get total gene population(genes in the module)
    geneFile=f"../coexpression_code/geneModule.txt"
    allgenesME=load_genes(geneFile)
    moduleGenes=allgenesME[module]
    
    #get GO class list
    gofile = "./genes/maize.B73.AGPv4.aggregate.gaf"
    golist=[]
    with open(gofile, "r") as gf:
        for line in gf:
            
            if "!" not in line[0]: #lines that start with ! are 1st info line and header
                golist.append(gpc.GOfile(line)) #use class from the created file
    
    #get v4 to v5 gene ID dictionary for translation
    translation_path="../data/annotation/genes_all.txt"
    translation=gpc.v4tov5_id(translation_path)
    
    #genes to GO terms to GO functions
    obo_file="./genes/go-basic.obo"
    obo_GOlist=GODag(obo_file)
    geneGov5=geneV5_GOterm(genelist, golist, translation)
    geneFunction=term_function(geneGov5, obo_GOlist)
    knG=list(geneFunction.keys()) #list  of known genes
    
    #Find unknown genes   
    ukw_genes=unknownGene(genelist, knG) #input is all of the comunity and all of the GO annotated genes
    print(f">Initial genes({len(genelist)})-GO term genes({len(geneGov5)})=unknown genes({len(ukw_genes)})\n") #missing genes appear in the counts tables but not in any reference
    print("Unknown gene analysis:\n")
    
    #Get genes closely correlated to the unknown ones
    corpath=f"./correlation_tables/{module}_geneCorrelation.txt"
    cor_df=pd.read_table(corpath, index_col=0)
    
    all_kcg={} #dictionary unknown gene:[known correlated genes, functions, GOterms]
    for g in genelist:

        if g in ukw_genes:#gene is unknown
            ug=g
            top_corGenes=closeToGene(ug, cor_df, n=n_corrGenes) #closely corr genes
            
            ug_kcg=[g for g in top_corGenes if g in knG] #extract known genes from the closely correlated genes list
            
            print(f"From {len(top_corGenes)} closely correlated genes to {ug}, {len(ug_kcg)} are known and belong to the community")
            
            if len(ug_kcg)==0: #If none of the correlated genes are known (or belong to the community)
                all_kcg[ug]=[ug_kcg, [], []] #determine empty dictionary
                continue
            
            #gene enrichment
            goEnrv5=GOEnrichmentStudy(moduleGenes, geneGov5, obo_GOlist, methods=['bonferroni', 'fdr_bh'], log=None) #v5 ID works as wll as V4 ID #remove log=None for detailed results
            
            results=goEnrv5.run_study(ug_kcg) #Run analysis on population(module) vs samples(known genes correlated to unknown)
            
            #get only significant go terms
            significant=[r for r in results if r.p_fdr_bh < 0.05] #use fdr_bh because better detail than provided by bonferroni
            #sort them by significance, get top n results
            sorted_sigGo=sorted(significant, key=lambda x: x.p_fdr_bh)[:n_sigGO]
            
            #view significant functions with the p-value
            '''for r in sorted_sigGo:
                print(f"{r.name}\t{r.p_fdr_bh:.2e}")'''
            
            
            #get correlated genes with signficant functions and GOterms for each unknown gene (genes woth no GO term or function presnet 2 empty lists as 1,2)
            all_kcg[ug]=[ug_kcg, [goTerm.name for goTerm in sorted_sigGo], [rec.GO for rec in sorted_sigGo]] #gene:[corr, function, goterm]
        
        elif g in knG: #gene is known

            all_kcg[g]=[["Known"], list(geneFunction[g]), list(geneGov5[g])] #gene:[[], function, goterm]
        
        print("#########################################################################################")
    
    
    print(f"Unknown genes={len(ukw_genes)}")
    return all_kcg
            

if __name__ == "__main__":
    genelist=['Zm00001eb199820', 'Zm00001eb088400', 'Zm00001eb391890', 'Zm00001eb216000', 'Zm00001eb390720', 'Zm00001eb350490', 'Zm00001eb279160', 'Zm00001eb416490', 'Zm00001eb264350', 'Zm00001eb070770', 'Zm00001eb305240', 'Zm00001eb228670', 'Zm00001eb104260', 'Zm00001eb176740', 'Zm00001eb373200', 'Zm00001eb114690', 'Zm00001eb408100', 'Zm00001eb154950', 'Zm00001eb049080', 'Zm00001eb361970', 'Zm00001eb318710', 'Zm00001eb305250', 'Zm00001eb400540', 'Zm00001eb054470', 'Zm00001eb156380', 'Zm00001eb114490', 'Zm00001eb169390', 'Zm00001eb074000', 'Zm00001eb100830', 'Zm00001eb080580', 'Zm00001eb168140', 'Zm00001eb402470', 'Zm00001eb190260', 'Zm00001eb366580', 'Zm00001eb376940', 'Zm00001eb057860', 'Zm00001eb376220', 'Zm00001eb233280', 'Zm00001eb204000', 'Zm00001eb127190', 'Zm00001eb233060', 'Zm00001eb137740', 'Zm00001eb328420', 'Zm00001eb325610', 'Zm00001eb315750', 'Zm00001eb082200', 'Zm00001eb086540', 'Zm00001eb079210', 'Zm00001eb350370', 'Zm00001eb114600', 'Zm00001eb082600', 'Zm00001eb348320', 'Zm00001eb420890', 'Zm00001eb156010', 'Zm00001eb376320', 'Zm00001eb098250', 'Zm00001eb221840', 'Zm00001eb341540', 'Zm00001eb128680', 'Zm00001eb392280', 'Zm00001eb210110', 'Zm00001eb071960', 'Zm00001eb065350', 'Zm00001eb226170', 'Zm00001eb088480', 'Zm00001eb030410', 'Zm00001eb283130', 'Zm00001eb299860', 'Zm00001eb074740', 'Zm00001eb056770', 'Zm00001eb015260', 'Zm00001eb251150', 'Zm00001eb123680', 'Zm00001eb076590', 'Zm00001eb093060', 'Zm00001eb320950', 'Zm00001eb043850', 'Zm00001eb110230', 'Zm00001eb088390', 'Zm00001eb211520', 'Zm00001eb150710', 'Zm00001eb017450', 'Zm00001eb079610', 'Zm00001eb114520', 'Zm00001eb231020', 'Zm00001eb269570', 'Zm00001eb258790', 'Zm00001eb376310', 'Zm00001eb140420', 'Zm00001eb376250', 'Zm00001eb264320', 'Zm00001eb351910', 'Zm00001eb231100', 'Zm00001eb224700', 'Zm00001eb289550', 'Zm00001eb433160', 'Zm00001eb376350', 'Zm00001eb404750', 'Zm00001eb090140', 'Zm00001eb179040', 'Zm00001eb002450', 'Zm00001eb170760', 'Zm00001eb379790', 'Zm00001eb253620', 'Zm00001eb201300', 'Zm00001eb090510', 'Zm00001eb001050', 'Zm00001eb275640', 'Zm00001eb035040', 'Zm00001eb124670', 'Zm00001eb208250', 'Zm00001eb393710', 'Zm00001eb345050', 'Zm00001eb364970', 'Zm00001eb243890', 'Zm00001eb030760', 'Zm00001eb421590', 'Zm00001eb376910', 'Zm00001eb297560', 'Zm00001eb267110', 'Zm00001eb000570', 'Zm00001eb118170', 'Zm00001eb192100', 'Zm00001eb292630', 'Zm00001eb022930', 'Zm00001eb004990', 'Zm00001eb402460', 'Zm00001eb035670', 'Zm00001eb422370', 'Zm00001eb191370', 'Zm00001eb104940', 'Zm00001eb411120', 'Zm00001eb084890', 'Zm00001eb313310', 'Zm00001eb159810', 'Zm00001eb196590', 'Zm00001eb028810', 'Zm00001eb360140', 'Zm00001eb376260', 'Zm00001eb140040', 'Zm00001eb020280', 'Zm00001eb231150', 'Zm00001eb229620', 'Zm00001eb044920', 'Zm00001eb173330', 'Zm00001eb117230', 'Zm00001eb356970', 'Zm00001eb207750', 'Zm00001eb184630', 'Zm00001eb082190', 'Zm00001eb404520', 'Zm00001eb082610', 'Zm00001eb046310', 'Zm00001eb402500', 'Zm00001eb230310', 'Zm00001eb433050', 'Zm00001eb254580', 'Zm00001eb195120', 'Zm00001eb114660', 'Zm00001eb082210', 'Zm00001eb345770', 'Zm00001eb201310', 'Zm00001eb257920', 'Zm00001eb169380', 'Zm00001eb425180', 'Zm00001eb137040', 'Zm00001eb405880', 'Zm00001eb231500', 'Zm00001eb012870', 'Zm00001eb376190', 'Zm00001eb070730']
    genelist=['Zm00001eb293680', 'Zm00001eb224870', 'Zm00001eb137020', 'Zm00001eb197030', 'Zm00001eb104850', 'Zm00001eb424940', 'Zm00001eb202550', 'Zm00001eb272420', 'Zm00001eb157270', 'Zm00001eb362170', 'Zm00001eb293670', 'Zm00001eb318200', 'Zm00001eb011640', 'Zm00001eb128830', 'Zm00001eb202560', 'Zm00001eb186790', 'Zm00001eb013250', 'Zm00001eb293660']
    #genelist=['Zm00001eb034780', 'Zm00001eb034790']
    module="violet"
    main(module, genelist, n_corrGenes=10, n_sigGO=10)
            