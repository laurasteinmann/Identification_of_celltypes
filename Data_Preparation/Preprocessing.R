# Save Root Atlas normaliced count data as csv-file
# Necessary Libraries
library(Seurat)
# Data Analysis
rc.integrated <- readRDS('/netscratch/dep_mercier/grp_novikova/laura/data/Single_Cell/GSE152766_Root_Atlas.rds') #Load data 
normalized_counts <- as.matrix(rc.integrated@assays$integrated@data) #extract normalized counts
celltypes_shahan <- as.matrix(rc.integrated@meta.data$celltype.anno) #extract annotation of Shahan et al.

write.csv(normalized_counts, "/netscratch/dep_mercier/grp_novikova/laura/data/Single_Cell/Root_Atlas_normalized_counts.csv", row.names=TRUE, col.names=TRUE) #Save counts as csv to import it with python 
write.csv(celltypes_shahan, "~/data/Single_Cell/Celltypes_Shahan.csv")