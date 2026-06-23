import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load matrix
matrix = pd.read_csv('gene_presence_absence_matrix.csv', index_col=0)

# Calculate number of genomes per gene
matrix['num_genomes'] = matrix.sum(axis=1)

# Sort genes by presence frequency
matrix_sorted = matrix.sort_values('num_genomes', ascending=False)

# Remove helper column for plotting
matrix_plot = matrix_sorted.drop(columns='num_genomes')

# Plot heatmap
plt.figure(figsize=(8, 12))
sns.heatmap(matrix_plot, cmap='viridis', yticklabels=False)
plt.xlabel('Genome')
plt.ylabel('Gene cluster')
plt.title('Gene Presence/Absence Across 4 Genomes')
plt.tight_layout()
plt.savefig('gene_presence_absence_heatmap.png', dpi=300)
plt.show()