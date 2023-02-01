## Taking the csv data save it as hdf5 and preprocess it for UMAP and HDBSCAN ##
## Import csv file with normalized counts and save as hdf5 file
import pandas as pd
data = pd.read_csv('/netscratch/dep_mercier/grp_novikova/laura/data/Single_Cell/Root_Atlas_normalized_counts.csv',index_col=0)
##### Get the corresponding celltypes
data = data.transpose()
data.to_hdf('/netscratch/dep_mercier/grp_novikova/laura/data/Single_Cell/RNA-Seq.h5', key='counts')

