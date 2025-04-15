#!/usr/bin/env python3

import pandas as pd
import numpy as np
from getSequence import getseq


genecode=pd.read_table("./geneModule.txt", comment="#",
                        sep = " ")
print(genecode.info())

print(getseq("Zm00001eb002950", uniprot_id=True, return_dict=True))

import os
import pandas as pd
import numpy as np
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.base import download_go_basic_obo
from goatools.associations import read_gaf

# Example function to parse and validate maize gene IDs
def validate_maize_ids(gene_ids):
    """Validate if the provided IDs follow Zea mays ID format."""
    valid_ids = []
    invalid_ids = []
    
    for gene_id in gene_ids:
        if gene_id.startswith('Zm') and len(gene_id) > 5:
            valid_ids.append(gene_id)
        else:
            invalid_ids.append(gene_id)
            
    print(f"Valid IDs: {len(valid_ids)}, Invalid IDs: {len(invalid_ids)}")
    return valid_ids, invalid_ids

# Example usage
maize_ids = ["Zm00001eb002950", "Zm00001eb034560", "Zm00001eb000123", "Invalid_ID"]
valid_ids, invalid_ids = validate_maize_ids(maize_ids)

##
import requests

def fetch_go_terms_for_gene(gene_id):
    """Fetch GO terms for a given gene ID from an API."""
    base_url = "https://maizegdb.org/api/v1/gene/"  # Example API endpoint
    response = requests.get(f"{base_url}{gene_id}")
    
    if response.status_code == 200:
        data = response.json()
        go_terms = data.get('go_terms', [])
        return go_terms
    else:
        print(f"Failed to retrieve data for {gene_id}: {response.status_code}")
        return []

##

# Download GO basic OBO file
obo_file = "go-basic.obo"
if not os.path.exists(obo_file):
    download_go_basic_obo(obo_file)

# Parse the GO OBO file
go = obo_parser.GODag(obo_file)

# Function to generate GO term summaries for genes
def summarize_go_terms(gene_ids, annotations, go_dag):
    """Generate a summary of GO terms for a list of genes."""
    results = {}
    
    for gene_id in gene_ids:
        if gene_id in annotations:
            go_terms = annotations[gene_id]
            go_info = []
            
            for go_id in go_terms:
                if go_id in go_dag:
                    term = go_dag[go_id]
                    go_info.append({
                        'go_id': go_id,
                        'name': term.name,
                        'namespace': term.namespace,
                        'depth': term.depth,
                        'is_obsolete': term.is_obsolete
                    })
            
            results[gene_id] = go_info
        else:
            results[gene_id] = []
            
    return results

##

import matplotlib.pyplot as plt
from collections import Counter

def plot_go_distribution(annotations, namespace_filter=None):
    """Plot the distribution of GO terms by namespace."""
    # Collect all GO namespaces
    namespaces = []
    for gene_id, go_terms in annotations.items():
        for term in go_terms:
            if namespace_filter is None or term['namespace'] == namespace_filter:
                namespaces.append(term['namespace'])
    
    # Count occurrences
    namespace_counts = Counter(namespaces)
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.bar(namespace_counts.keys(), namespace_counts.values())
    plt.xlabel('GO Namespace')
    plt.ylabel('Count')
    plt.title('Distribution of GO Terms by Namespace')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

##