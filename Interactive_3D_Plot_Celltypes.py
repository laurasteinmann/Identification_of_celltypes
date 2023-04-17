import os
import h5py
import pandas as pd
import numpy as np
import plotly.express as px
# Plot 3D Cluster Classification on UMAP
os.chdir('../../Single_Cell_Data/')
filename = 'UMAP/umap3_24.h5'
with h5py.File(filename, "r") as f:
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    data = pd.DataFrame(f[a_group_key])

data.columns = ["UMAP1", "UMAP2", "UMAP3"]
cluster = pd.read_csv("HDBSCAN_Best/Labels_UMAP4_20_99_79.csv", header=None)
data["Cluster"] = np.array(cluster.iloc[:, 0], dtype=int)
fig = px.scatter_3d(data, x='UMAP1', y='UMAP2', z='UMAP3', color='Cluster')
fig.write_html("3D-Plots/3D_UMAP_Bestrun_Visualisation.html")
fig.show()