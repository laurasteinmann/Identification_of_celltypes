from Functions import sort_matrix, normalizing_matrix
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
# needed for importing the functions from Functions don't know exactly why but solves import problems
sys.path.insert(1, '')
# Creating comparative tables of biological and machine clustering
os.chdir('../../Single_Cell_Data/')
results = pd.read_csv("Parameters_Bestratio_old.csv", index_col=0)
results.shape
results
runs = range(0, results.shape[0])
for run in runs:
    parameters = results.to_dict('records')[run]
    print(parameters)
    labels_file = "HDBSCAN_Best/Labels_UMAP" + str(int(parameters["UMAP_Dim"])) + "_" + str(int(parameters['UMAP_n_neighbors'])) + "_" + str(
        int(parameters["HDBSCAN_min_samples"])) + "_" + str(int(parameters["HDBSCAN_cluster_size"])) + ".csv"
    output_file = "Comparison_Tables/Comparison_" + str(int(parameters["UMAP_Dim"])) + "_" + str(int(parameters['UMAP_n_neighbors'])) + "_" + str(
        int(parameters["HDBSCAN_min_samples"])) + "_" + str(int(parameters["HDBSCAN_cluster_size"])) + ".csv"
    cluster = np.loadtxt(labels_file, delimiter=",")
    celltype_labels = pd.read_csv('Celltypes_Shahan.csv', index_col=0)
    celltype_labels = np.array(celltype_labels['V1'])
    celltypes = np.unique(celltype_labels)
    comparison = pd.DataFrame(0, columns=np.unique(cluster), index=celltypes)
    comparison.columns = comparison.columns.astype(int)
    if celltype_labels.shape[0] == cluster.shape[0]:
        for i in np.arange(0, len(cluster)):
            column = int(cluster[i])
            row = celltype_labels[i]
            comparison.loc[row, column] = comparison.loc[row, column] + 1
    else:
        print("Problem occurred during clustering. There are different numbers of cells")

    comparison = comparison.transpose()
    sorted_comparison = sort_matrix(comparison)
    normalized_matrix = normalizing_matrix(sorted_comparison)
    normalized_matrix = normalized_matrix.transpose()
    normalized_matrix.to_csv(output_file)

# Plot comparison table as heatmap
results = pd.read_csv("Parameters_Bestratio_old.csv", index_col=0)
results.shape
runs = range(0, results.shape[0])
for run in runs:
    parameters = results.to_dict('records')[run]
    comparison_file = "Comparison_Tables/Comparison_" + str(int(parameters["UMAP_Dim"])) + "_" + str(int(
        parameters['UMAP_n_neighbors'])) + "_" + str(int(parameters["HDBSCAN_min_samples"])) + "_" + str(int(parameters["HDBSCAN_cluster_size"])) + ".csv"
    comparison = pd.read_csv(comparison_file, index_col=0)
    output_file = "Heatmaps/Results_Comparative_Heatmap_" + str(int(parameters["UMAP_Dim"])) + "_" + str(int(parameters['UMAP_n_neighbors'])) + "_" + str(
        int(parameters["HDBSCAN_min_samples"])) + "_" + str(int(parameters["HDBSCAN_cluster_size"])) + ".png"
    fig = plt.subplots(1, 1, figsize=(12, 8))
    heat_map = sns.heatmap(comparison, cmap='inferno', cbar_kws={
                           'label': 'Percentage of cells of one celltype'})
    heat_map.set_xticklabels(heat_map.get_xticklabels(), rotation=0)
    heat_map.set_yticklabels(heat_map.get_yticklabels(), rotation=35)
    plt.title('Best Randindex:' + str(parameters['Rand_Index'])
              + "__" + "Clusters:" + str(int(parameters['Clusters'])))
    plt.xlabel('HDBSCAN Clusters')
    plt.ylabel('Annotated Cell types')
    plt.tight_layout()
    plt.savefig(output_file, dpi=700, facecolor='w')
    plt.close()