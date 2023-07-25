

library(Seurat)
combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')

# data_counts <- as.matrix(combined_obj@assays$integrated@counts)
# rank.range <- seq(5, 20)
# for (rank in rank.range){
# 	print(rank)
# 	res <- NMF::nmf(data_counts, rank)
# 	saveRDS(res, paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Integrated_Rank_',rank,'.rds'))
# }


data_counts <- as.matrix(combined_obj@assays$integrated@data)
data_counts[data_counts<0] <- 0
rank.range <- seq(5, 20)
for (rank in rank.range){
	print(rank)
	res <- NMF::nmf(data_counts, rank)
	saveRDS(res, paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
}


# NMF::summary(res)
# a <- NMF::predict(res)


