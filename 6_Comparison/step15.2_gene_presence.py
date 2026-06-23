import pandas as pd

# Load cluster table
clusters = pd.read_csv('clusters.tsv', sep='\t', header=None, names=['rep','member'])

# Get all genes in all genomes
all_genes = pd.concat([clusters['rep'], clusters['member']])

# Extract unique genome prefixes
genomes = list(set([g.split('_')[0] for g in all_genes]))
genomes.sort()  # optional, for consistent column order

# Get all clusters (unique representatives)
cluster_ids = clusters['rep'].unique()

# Initialize matrix
matrix = pd.DataFrame(0, index=cluster_ids, columns=genomes)

# Fill matrix
for _, row in clusters.iterrows():
    for gene in [row['rep'], row['member']]:
        genome = gene.split('_')[0]
        matrix.loc[row['rep'], genome] = 1

# Save
matrix.to_csv('gene_presence_absence_matrix.csv')