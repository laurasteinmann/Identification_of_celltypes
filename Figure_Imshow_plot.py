# Figure I Supplement Imshow-Plot

## Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## Set Working directory
os.chdir('../../Single_Cell_Data/')
results = pd.read_csv("Parameters_old.csv", index_col=0)
results_dim4 = results[results["UMAP_Dim"]==4]
output_file = "Imshow_Plot.png"

hdbscan_samples = np.arange(3, 101)
data_list = []
n_rows = 4
umap_n_neighbors = results_dim4['UMAP_n_neighbors']
umaps = list(np.unique(umap_n_neighbors))
umaps = umaps[1:]
n_columns = int(np.ceil(len(umaps)/n_rows))
x_label = "Cluster Size"
y_label = "Min Samples"

fig, ax = plt.subplots(n_rows, n_columns, sharex=True,
                       sharey=True, figsize=(12, 8))
for plot_idx in range(len(umaps)):
    umap = umaps[plot_idx]
    umap2 = results_dim4[results_dim4['UMAP_n_neighbors'] == umap]
    umap2_full = umap2.iloc[:, 2:]
    data = pd.DataFrame(0, index=hdbscan_samples, columns=hdbscan_samples)
    for sample in hdbscan_samples:
        short = umap2_full[umap2_full['HDBSCAN_min_samples'] == sample]
        data[sample] = np.asarray(short['Rand_Index'])
    data = data.sort_index(axis=0, ascending=False)
    im = ax.flatten()[plot_idx].imshow(data, cmap='inferno',
                                       vmin=0, vmax=1, extent=[3, 99, 3, 99])
    title = 'UMAP ' + str(int(umap))
    ax.flatten()[plot_idx].set_title(title)
    
    
for row in range(0, n_rows):
    ax[row,0].set_ylabel(y_label)
    
for column in range(0, n_columns):
    ax[3,column].set_xlabel(x_label)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.legend()
plt.savefig(output_file, dpi=700, facecolor='w')
plt.show()