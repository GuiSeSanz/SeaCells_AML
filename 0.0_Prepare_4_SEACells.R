library(Seurat)
library(SeuratData)
library(SeuratDisk)


DATA_PATH  <-  '/home/sevastopol/data/gserranos/SEACells_AML/Data'
for (file in list.files(DATA_PATH, pattern = "\\.rds$", full.names = TRUE)){
	if (grepl("h5Seurat", file)){
		print(paste0("Skipping ", file))
		next
	}else{
		print(paste0("Processing ", file))
		data_aml <- readRDS(file)
		new_file <- stringr::str_replace(file, '\\.rds$', '.h5Seurat')
		data_aml[['SCT']] <- NULL
		DefaultAssay(data_aml) <- "RNA"
		# data_aml <- NormalizeData(data_aml)
		# data_aml <- ScaleData(data_aml, features=rownames(data_aml))
		# data_aml <- FindVariableFeatures(object = data_aml)
		# data_aml <- RunPCA(data_aml)
		# data_aml <- RunUMAP(data_aml, dims=1:50)
		SaveH5Seurat(data_aml, filename = new_file, overwrite=TRUE)
		Convert(new_file, dest = "h5ad", overwrite=TRUE)
	}
}





