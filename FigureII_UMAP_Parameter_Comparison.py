# Script for generating Figure II 
## Libraries
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
## Set Working directory
os.chdir('../../Single_Cell_Data/')
## Data
results = pd.read_csv("Parameters_old.csv", index_col=0)
results_dim4 = results[results["UMAP_Dim"]==4]
results_dim4_n24 = results_dim4[results_dim4["UMAP_n_neighbors"]==24]

x_labels_umap_dim = np.array(np.unique(results["UMAP_Dim"]), dtype=int)
x_labels_neighbors = np.array(np.unique(results["UMAP_n_neighbors"]), dtype=int)

fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 3, 1]}, sharey=True, figsize=(12, 4))
sns.violinplot(ax = ax[0], data=results, x = "UMAP_Dim", y = "Rand_Index")
ax[0].set_xticklabels(x_labels_umap_dim)
ax[0].set_ylabel("Rand Index")
sns.boxplot(ax = ax[1], data=results, x = "UMAP_n_neighbors", y = "Rand_Index")
ax[1].set_xticklabels(x_labels_neighbors)
sns.scatterplot(ax = ax[2], data=results, x = "Clusters", y = "Rand_Index", color="C0", linewidth=0.1, s = 10)
output_file = "FigureII.png"
plt.tight_layout()
plt.savefig(output_file, dpi=700, facecolor='w')
plt.show()
plt.close()