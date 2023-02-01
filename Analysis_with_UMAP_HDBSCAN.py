####### Use Computational Methods to classify clusters which represents #######
##### the biological celltypes and validate them with statistical results #####

######################### Libraries and Functions #############################
import sys
import pandas as pd
import seaborn as sns
import umap
import h5py
import hdbscan
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score

############################ Load Data #########################################
data = pd.read_hdf('/netscratch/dep_mercier/grp_novikova/laura/data/Single_Cell/RNA-Seq.h5', key='counts')
celltype_labels = pd.read_csv('Celltypes_Shahan.csv', index_col=0)
################### Prepare Celltype Labels for Comparison #####################
celltype_labels.insert(1, "Number",celltype_labels.iloc[:,0])
celltype_labels.columns = ['Celltype', 'Number']
names = celltype_labels.Celltype.unique()
for i in np.arange(0, len(names)):
    celltype_labels['Number'] = celltype_labels['Number'].replace([names[i]], i+1)

celltype_values = celltype_labels.loc[:,'Number']
control_clusters = celltype_values.values
###### Prepare Pandas Dataframe for saving Parameter settings and Results ######
columns = ["UMAP_Dim", "UMAP_n_neighbors", "HDBSCAN_min_samples",
           "HDBSCAN_cluster_size", "Clusters", "Rand_Index"]
results = pd.DataFrame(columns=columns)
############################# Run UMAP #########################################
umap_dim = 2
umap_neighbors = 2 #between 2 and 13 can be up to 200
umap_embedding = umap.UMAP(
        n_neighbors = umap_neighbors,
        min_dist = 0.0,
        n_components = umap_dim,
        random_state = 42,
        ).fit_transform(data)
        #Save UMAP in hdf5
umap_file = "UMAP/umap"+str(umap_dim)+"_"+str(umap_neighbors)+".h5"
h5f = h5py.File(umap_file,'w')
h5f.create_dataset('umap', data = umap_embedding)
h5f.close()
######################### HDBSCAN parameter ####################################
for j in np.arange(3, 101):
    parameter_samples = j
    for k in np.arange(3, 101):
        parameter_cluster = k
        labels = hdbscan.HDBSCAN(
                min_samples=int(parameter_samples),
                min_cluster_size=int(parameter_cluster),
                ).fit_predict(umap_embedding)
        clustered = (labels >= 0)
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
###################### Calculate Rand-Index ####################################
        rand_index = adjusted_rand_score(control_clusters, labels)
################## Add new scores to Results table #############################
        new_row = {"UMAP_Dim":umap_dim, "UMAP_n_neighbors":umap_neighbors,
                        "HDBSCAN_min_samples":parameter_samples,
                        "HDBSCAN_cluster_size": parameter_cluster,
                        "Clusters": n_clusters_, "Rand_Index": rand_index}
        results = results.append(new_row, ignore_index=True)
        results_file = "Parameters" + str(umap_dim) + "_" + str(umap_neighbors) + ".csv"
        results.to_csv(results_file)
        labels_file = "HDBSCAN/Labels_UMAP" + str(umap_dim) + "_" + str(umap_neighbors) + "_" + str(parameter_samples) + "_" + str(parameter_cluster) + ".csv"
            # fig, (ax1, ax2) = plt.subplots(1, 2)
            # fig.suptitle('UMAP3 '+'HDBSCAN '+ 'Samples '+str(parameter_samples)+'Clusters'+ str(parameter_cluster)+'Number of Clusters:'+str(n_clusters_))
            # ax1.scatter(umap_embedding[clustered, 0],umap_embedding[clustered, 1],
            #             c=labels[clustered],s=0.2,cmap='turbo')
            # ax2.scatter(umap_embedding[clustered, 0],
            #             umap_embedding[clustered, 1],
            #             c=labels[clustered],
            #             s=0.2,
            #             cmap='turbo')
            # ax1.set_ylabel("UMAP_2")
            # ax1.set_xlabel("UMAP1")
            # ax1.set_ylabel("UMAP_2")
            # ax2.set_xlabel("UMAP_3")
            # output_file = "HDBSCAN/UMAP" + str(umap_dim) +"_"+ str(umap_neighbors) +"_HDBSCAN_"+ str(parameter_samples)+"_"+str(parameter_cluster)+".png"
            # plt.savefig(output_file, dpi=700)
            #plt.close()
        np.savetxt(labels_file, labels, delimiter=",")