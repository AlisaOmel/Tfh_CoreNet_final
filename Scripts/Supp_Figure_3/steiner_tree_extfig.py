import pandas as pd
import numpy as np
import csv
import h5py
import pickle
import random
random.seed(42)
import ast
import networkx as nx
from networkx.algorithms.approximation import steiner_tree
from tqdm import tqdm
from itertools import combinations

ppi = h5py.File('/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5', 'r')
hint_data = ppi['edges']

#with open("/ix/djishnu/Alisa/Tfh/Vinuesa_updated_list.txt", "r") as file:
#    Vinuesa_review = file.read().strip()  # Read and remove any leading/trailing spaces
    
# Convert string to set
#Vinuesa_review_gene_set = ast.literal_eval(Vinuesa_review)

vinuesa_df = pd.read_csv("../../Data/literature_sets/vinuesa_list.csv")
Vinuesa_review_gene_set = set(vinuesa_df['Genes'].dropna().astype(str).str.strip())

### Make a dataframe with HINT interactions and the Gene IDs
def make_hint_df(hint_data):
    gene_sets = []
    for entry in hint_data:
        prot_1 = entry[0].decode('utf-8')
        prot_2 = entry[1].decode('utf-8')
        
        gene_sets.append([prot_1, prot_2])
    return gene_sets

hint = make_hint_df(hint_data)
hint_df = pd.DataFrame(hint)

G= nx.from_pandas_edgelist(hint_df, 0,1)

terminals = [g for g in Vinuesa_review_gene_set if g in G]
missing = Vinuesa_review_gene_set - set(terminals)
if missing:
    print(f"⚠️ Missing {len(missing)} genes from HINT: {missing}")


steiner_subgraph = steiner_tree(G, terminals)

# 6. Export or visualize
nx.write_edgelist(steiner_subgraph, "steiner_subnetwork.edgelist")
print(f"Steiner network: {steiner_subgraph.number_of_nodes()} nodes, {steiner_subgraph.number_of_edges()} edges")


# === 1. Export to Cytoscape .sif format ===
with open("vinuesa_steiner_network.sif", "w") as f:
    for u, v in steiner_subgraph.edges():
        f.write(f"{u}\tpp\t{v}\n")

# === 2. Create and export node attribute file ===
nodes = list(steiner_subgraph.nodes())
labels = pd.DataFrame({
    "node_id": nodes,
    "in_gene_set": ["TRUE" if n in Vinuesa_review_gene_set else "FALSE" for n in nodes]
})
labels.to_csv("vinuesa_steiner_node_attributes.txt", sep="\t", index=False)

print("Files written:\n- steiner_network.sif\n- steiner_node_attributes.txt")
