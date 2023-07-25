

library(Seurat)
library(ggplot2)
library(cluster)  # Required for silhouette function
library(factoextra)  # Required for fviz_nbclust function

spectral_clustering <- function(X, # matrix of data points
																nn = 10, # the k nearest neighbors to consider
																n_eig = 2) # m number of eignenvectors to keep
	{
	mutual_knn_graph <- function(X, nn = 10)
	{
		D <- as.matrix( dist(X) ) # matrix of euclidean distances between data points in X
		
		# intialize the knn matrix
		knn_mat <- matrix(0,
											nrow = nrow(X),
											ncol = nrow(X))
		
		# find the 10 nearest neighbors for each point
		for (i in 1: nrow(X)) {
			neighbor_index <- order(D[i,])[2:(nn + 1)]
			knn_mat[i,][neighbor_index] <- 1 
		}

		# Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1 
		knn_mat <- knn_mat + t(knn_mat) # find mutual knn
		
		knn_mat[ knn_mat == 2 ] = 1
		
		return(knn_mat)
	}

	graph_laplacian <- function(W, normalized = TRUE)
	{
		stopifnot(nrow(W) == ncol(W)) 
		
		g = colSums(W) # degrees of vertices
		n = nrow(W)
		
		if(normalized)
		{
			D_half = diag(1 / sqrt(g) )
			return( diag(n) - D_half %*% W %*% D_half )
		}
		else
		{
			return( diag(g) - W )
		}
	}
	print(paste0('Calculating similarity matrix for a matrix of size: ', c(dim(X))))
	W = mutual_knn_graph(X) # 1. matrix of similarities
	print(paste0('Calculating graph laplacian for a matrix of size: ', dim(W)))
	L = graph_laplacian(W) # 2. compute graph laplacian
	print(paste0('Calculating eigenvalues for a matrix of size: ', dim(L)))
	ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
	n = nrow(L)
	print(paste0('Calculating eigenvectors for a matrix of size: ', dim(ei$vectors)))
	return(ei$vectors[,(n - n_eig):(n - 1)]) # return the eigenvectors of the n_eig smallest eigenvalues

}

create_folder <- function(folder_path) {
	if (!file.exists(folder_path)) {
		dir.create(folder_path, recursive = TRUE)
		cat("Folder created:", folder_path, "\n")
	} else {
		cat("Folder already exists:", folder_path, "\n")
	}
	return(folder_path)
}

calculate_sparsity <- function(counts){
	sparsity <- 1 - sum(counts > 0) / (dim(counts)[1] * dim(counts)[2])
	return(sparsity)
}


get_kmeans_results <- function(data, k){
	clusters <- kmeans(X_sc, k)$cluster
	sil_score <- silhouette(clusters, dist(data))[,'sil_width']
	return(list(clusters=clusters, sil_score=sil_score))
}

DATA_FOLDER = '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells'
PLOT_PATH = '/home/sevastopol/data/gserranos/SEACells_AML/Plots/Single_cell_VS_SEACells'
AML_SAMPLES = grep('aml', list.dirs(DATA_FOLDER, full.names=FALSE), value=TRUE)


all_metacells_objects <- NULL
Sparsity_results <- data.frame(Sample=NULL, Sparsity_sc=NULL, Sparsity_SEACells=NULL)
for (aml_sample in AML_SAMPLES){
	print(aml_sample)
	Metacells          <- read.csv(  paste0(DATA_FOLDER, '/', aml_sample, '/AML_seacells_data_colnames.csv'), sep='\t', header=TRUE, row.names=1)
	Metacells_metadata <- read.table(paste0(DATA_FOLDER, '/', aml_sample, '/AML_seacells_metadata.csv'), sep=' ', header=TRUE)
	Metacells_coords   <- read.table(paste0(DATA_FOLDER, '/', aml_sample, '/AML_seacells_coords_colnames.csv'), sep='\t', header=TRUE, row.names=1)
	# leiden             <- read.table(paste0(DATA_FOLDER, '/', aml_sample, '/AML_seacells_LeidenClustering.csv'),  sep='\t', header=TRUE)
	metacells_Seurat <- CreateSeuratObject(counts = as.matrix(t(Metacells)), project = aml_sample, min.cells = 3, min.features = 200)
	rawCounts <- metacells_Seurat@assays$RNA@counts
	counts_SEACells <- t(as.matrix(read.csv(paste0(DATA_FOLDER, '/', aml_sample, '/AML_seacells_RAW_data_colnames.csv'), sep='\t', header=TRUE, row.names=1)))
	metacells_Seurat <- SeuratObject::SetAssayData(metacells_Seurat, slot='counts', new.data = counts_SEACells[rownames(metacells_Seurat), ], assay = "RNA")
	rownames(Metacells_coords) <- Metacells_coords$SEACell
	Metacells_coords <- Metacells_coords[,c('X0', 'X1')]
	colnames(Metacells_coords) <- paste0("SEACell_", 1:2)
	Metacells_coords <- as.matrix(Metacells_coords)
	metacells_Seurat[["umap_metacells"]] <- CreateDimReducObject(embeddings = Metacells_coords, key = "SEACell_", assay = DefaultAssay(metacells_Seurat))
	aml_sc <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/seurat_', aml_sample, '.rds'))
	origin_coordinates <-  read.csv(paste0(DATA_FOLDER, '/', aml_sample, '/original_umap_coords.csv'), sep='\t', header=TRUE, row.names=1)
	aml_sc[["umap_origin"]] <-  CreateDimReducObject(embeddings =  as.matrix(origin_coordinates), key = "originumap_", assay = DefaultAssay(aml_sc))
	plot_path_tmp <- create_folder(paste0(PLOT_PATH, '/', aml_sample))
	coords_SEACell <- FetchData(metacells_Seurat, c('SEACell_1', 'SEACell_2'))
	coords_sc      <- FetchData(aml_sc, c('originumap_0', 'originumap_1'))
	
	# Comparison of the different umap coordinates
	pdf(paste0(plot_path_tmp, '/Umap_Comparison.pdf'))
		print(cowplot::plot_grid(ggplot(coords_SEACell, aes(x=coords_SEACell[,1], y=coords_SEACell[,2])) + geom_point() + theme_classic(),
		ggplot(coords_sc,      aes(x=coords_sc[,1], y=coords_sc[,2])) + geom_point() + theme_classic(), ncol=2))
		print(ggplot() +
		geom_point(data= coords_sc, mapping=aes(x=coords_sc[,1], y=coords_sc[,2]), color= '#AEF098', size=0.3) +
		geom_point(data= coords_SEACell, mapping=aes(x=coords_SEACell[,1], y=coords_SEACell[,2]), color='#962493', size=0.6) + 
		theme_classic())
	dev.off()

	# Sparsity check
	SEACell_counts <- as.matrix(metacells_Seurat@assays$RNA@counts)
	sc_counts <- as.matrix(aml_sc@assays$RNA@counts)
	Sparsity_results <- rbind(Sparsity_results, data.frame(Sample=aml_sample, 
												Sparsity_sc=calculate_sparsity(sc_counts),
												 Sparsity_SEACells=calculate_sparsity(SEACell_counts)))
	metacells_Seurat <- RenameCells(metacells_Seurat, new.names = paste0(aml_sample, '_', colnames(SEACell_counts)))
	if(is.null(all_metacells_objects)){
		all_metacells_objects <- metacells_Seurat
	}else{
		all_metacells_objects <- merge(all_metacells_objects, metacells_Seurat)
	}
}


palette <- c(
		"#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
		"#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
		"#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
		"#f3e1eb", "#f6c4e1", "#f79cd4", "#7f7f7f", "#c7c7c7", "#1CE6FF", "#336600"
)
pdf("/home/sevastopol/data/gserranos/SEACells_AML/Plots/Sparsity_check.pdf")
plotter <- setNames(reshape2::melt(Sparsity_results), c('Sample', 'Assay', 'Sparsity'))
ggplot(plotter, aes(x=Assay, y=Sparsity, group=Sample)) +  geom_line() + geom_point(aes(color=Sample), size=5, alpha=0.8) + theme_classic() + 
scale_color_manual(values=palette) 
dev.off()



# apply to the merged object
# saveRDS(all_metacells_objects, '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells.rds')
all_metacells_objects <- readRDS('/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells.rds')
all_metacells_objects <- FindVariableFeatures(all_metacells_objects, selection.method = "vst", nfeatures = 3000)

metacells_normalized_data <- as.matrix(all_metacells_objects@assays$RNA@data)[c(all_metacells_objects@assays$RNA@var.features),]
X_sc <- spectral_clustering(t(metacells_normalized_data))
# run kmeans on the 2 eigenvectors
X_sc_kmeans <- kmeans(X_sc, 3)

k_range <- 2:15 
# Calculate silhouette scores for different K values
silhouette_scores <- sapply(k_range, function(k) {
  kmeans_res <- kmeans(X_sc, centers = k, nstart = 10)  # Perform K-means clustering
  silhouette_avg <- mean(silhouette(kmeans_res$cluster, dist(X_sc))[,'sil_width'] ) # Calculate average silhouette width
  return(silhouette_avg)
})

# Plot the silhouette scores
pdf(paste0(PLOT_PATH, '/Silhouette_analysis_Spectral.pdf'))
plot(k_range, silhouette_scores, type = "b", xlab = "Number of Clusters (K)", ylab = "Silhouette Score")
fviz_nbclust(X_sc, kmeans, method = "silhouette")
dev.off()

all_metacells_objects <- ScaleData(all_metacells_objects, verbose = FALSE)
all_metacells_objects <- RunPCA(all_metacells_objects, features = VariableFeatures(object = all_metacells_objects), verbose = FALSE)
all_metacells_objects <- FindNeighbors(all_metacells_objects, dims = 1:30)
all_metacells_objects <- RunUMAP(all_metacells_objects, reduction = "pca", dims = 1:30)



coords <- FetchData(all_metacells_objects, c('UMAP_1', 'UMAP_2', 'Sample'))

tmp <- get_kmeans_results(X_sc, 3)
coords$k_3 <- factor(tmp$clusters)
coords$silhouette_K3 <-  tmp$sil_score

tmp <- get_kmeans_results(X_sc, 8)
coords$k_8 <- factor(tmp$clusters)
coords$silhouette_K8 <-  tmp$sil_score

pdf(paste0(PLOT_PATH, '/UMAP_merged_SpectralClustering_performance.pdf'))
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=Sample)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
scale_color_manual(values=palette) + theme(legend.position="bottom")

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral") + theme(legend.position="bottom"), 
ncol=2)

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral")+ theme(legend.position="bottom"), 
ncol=2)
dev.off()

####################################################################################
# RPCA
####################################################################################

all_metacells_objects$Sample <- all_metacells_objects$orig.ident
all_metacells_list <- SplitObject(all_metacells_objects, split.by='Sample')

all_metacells_list <- lapply(X = all_metacells_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = all_metacells_list)
all_metacells_list <- lapply(X = all_metacells_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = all_metacells_list, anchor.features = features, reduction = "rpca", k.filter=50)
combined_obj <- IntegrateData(anchorset = anchors, k.weight=25)
DefaultAssay(combined_obj) <- "integrated"


combined_obj <- ScaleData(combined_obj, verbose = FALSE)
combined_obj <- RunPCA(combined_obj, npcs = 30, verbose = FALSE)
combined_obj <- RunUMAP(combined_obj, reduction = "pca", dims = 1:30)
combined_obj <- FindNeighbors(combined_obj, reduction = "pca", dims = 1:30)
combined_obj <- FindClusters(combined_obj, resolution = 0.5)

plotter <- FetchData(combined_obj, c('UMAP_1', 'UMAP_2', 'Sample', 'integrated_snn_res.0.5'))
pdf(paste0(PLOT_PATH, '/UMAP_combined.pdf'))
cowplot::plot_grid(
ggplot(plotter, aes(x=UMAP_1, y=UMAP_2, color=Sample)) + geom_point(size=0.5, alpha=0.8) + 
theme_classic() + scale_color_manual(values=palette) + theme(legend.position="bottom"),
ggplot(plotter, aes(x=UMAP_1, y=UMAP_2, color=integrated_snn_res.0.5)) + geom_point(size=0.5, alpha=0.8) + 
theme_classic()+ scale_color_manual(values=palette) + theme(legend.position="bottom"),
ncol=2)
dev.off()

# saveRDS(combined_obj, '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated.rds')
combined_obj <- readRDS('/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated.rds')
##### Apply the Spectral Clustering calcultaed on the merged to theintegrated data

DefaultAssay(combined_obj) <- 'RNA'
coords <- FetchData(combined_obj, c('UMAP_1', 'UMAP_2', 'Sample'))


metacells_normalized_data <- as.matrix(all_metacells_objects@assays$RNA@data)[c(all_metacells_objects@assays$RNA@var.features),]
X_sc <- spectral_clustering(t(metacells_normalized_data))

tmp <- get_kmeans_results(X_sc, 3)
coords$k_3 <- factor(tmp$clusters)
coords$silhouette_K3 <-  tmp$sil_score

tmp <- get_kmeans_results(X_sc, 8)
coords$k_8 <- factor(tmp$clusters)
coords$silhouette_K8 <-  tmp$sil_score

pdf(paste0(PLOT_PATH, '/UMAP_RPCA_SpectralClustering_performance.pdf'))
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=Sample)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
scale_color_manual(values=palette) + theme(legend.position="bottom")

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral") + theme(legend.position="bottom"), 
ncol=2)

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral")+ theme(legend.position="bottom"), 
ncol=2)
dev.off()






# TOO MUCH BATCH EFFECT!

DefaultAssay(combined_obj) <- 'integrated'
metacells_integrated_data <- as.matrix(combined_obj@assays$integrated@data)[c(combined_obj@assays$integrated@var.features),]
X_sc <- spectral_clustering(t(metacells_integrated_data), nn=10, n_eig=2)
# run kmeans on the 2 eigenvectors

k_range <- 2:15 
# Calculate silhouette scores for different K values
silhouette_scores <- sapply(k_range, function(k) {
  kmeans_res <- kmeans(X_sc, centers = k, nstart = 10)  # Perform K-means clustering
  silhouette_avg <- mean(silhouette(kmeans_res$cluster, dist(X_sc))[,'sil_width'] ) # Calculate average silhouette width
  return(silhouette_avg)
})

# Plot the silhouette scores
pdf(paste0(PLOT_PATH, '/Silhouette_analysis_Spectral_Integrated.pdf'))
plot(k_range, silhouette_scores, type = "b", xlab = "Number of Clusters (K)", ylab = "Silhouette Score")
fviz_nbclust(X_sc, kmeans, method = "silhouette")
dev.off()

coords <- FetchData(combined_obj, c('UMAP_1', 'UMAP_2', 'Sample'))
tmp <- get_kmeans_results(X_sc, 3)
coords$k_3 <- factor(tmp$clusters)
coords$silhouette_K3 <-  tmp$sil_score

tmp <- get_kmeans_results(X_sc, 8)
coords$k_8 <- factor(tmp$clusters)
coords$silhouette_K8 <-  tmp$sil_score

pdf(paste0(PLOT_PATH, '/UMAP_RPCA_SpectralClustering_performance_Integrated.pdf'))
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=Sample)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
scale_color_manual(values=palette) + theme(legend.position="bottom")

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K3)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral") + theme(legend.position="bottom"), 
ncol=2)

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_K8)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_distiller(palette = "Spectral")+ theme(legend.position="bottom"), 
ncol=2)
dev.off()



pdf(paste0(PLOT_PATH, '/UMAP_RPCA_SpectralClusteringVslouvain_Integrated.pdf'))
DefaultAssay(combined_obj) <- 'integrated'
coords <- FetchData(combined_obj, c('UMAP_1', 'UMAP_2', 'Sample', 'integrated_snn_res.0.5'))
tmp <- get_kmeans_results(X_sc, 11)
coords$k_11 <- factor(tmp$clusters)
coords$silhouette_K11 <-  tmp$sil_score

cowplot::plot_grid(
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_11)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette) + theme(legend.position="bottom"),
	ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=integrated_snn_res.0.5)) + geom_point(size=0.5, alpha=0.8) + theme_classic() +
	scale_color_manual(values=palette)+ theme(legend.position="bottom"), 
ncol=2)

dev.off()



# silhouette_avg <- silhouette(tmp$clusters, dist(X_sc))[,'sil_width'] 

DefaultAssay(combined_obj) <- 'integrated'
coords <- FetchData(combined_obj, c('UMAP_1', 'UMAP_2', 'Sample', 'integrated_snn_res.0.5'))
# coords$silhouette_spectral <- silhouette_avg
DIST_METACELLS_INTEGRATED_DATA <- dist(t(metacells_integrated_data))
coords$silhouette_louvain_data <- silhouette(as.numeric(coords$integrated_snn_res.0.5), DIST_METACELLS_INTEGRATED_DATA)[,'sil_width']

# based on the eigenvectors of the graph laplacian

pdf(paste0(PLOT_PATH, '/UMAP_RPCA_SpectralClusteringVslouvain_Integrated_Silhouette.pdf'))
for (nn in 1:11){
	print(nn)

	X_sc <- spectral_clustering(t(metacells_integrated_data), nn=nn, n_eig=11)
	tmp <- get_kmeans_results(X_sc, 11)
	coords$k_11 <- factor(tmp$clusters)
	coords$silhouette_Spectral_eigen <-  tmp$sil_score
	coords$silhouette_Spectral_data <- silhouette(tmp$clusters,DIST_METACELLS_INTEGRATED_DATA)[,'sil_width']


	coords$dunn_louvain_data   <- clValid::dunn(distance = DIST_METACELLS_INTEGRATED_DATA, as.numeric(coords$integrated_snn_res.0.5), Data = NULL, method = "euclidean")
	coords$dunn_spectral_eigen <- clValid::dunn(distance = dist(X_sc), tmp$clusters, Data = NULL, method = "euclidean")
	coords$dunn_spectral_data  <- clValid::dunn(distance = DIST_METACELLS_INTEGRATED_DATA, tmp$clusters, Data = NULL, method = "euclidean")

	coords$davies_bouldin_louvain_data   <- clusterSim::index.DB(t(metacells_integrated_data), as.numeric(coords$integrated_snn_res.0.5), DIST_METACELLS_INTEGRATED_DATA, centrotypes="centroids")$DB
	coords$davies_bouldin_spectral_eigen <- clusterSim::index.DB(X_sc, tmp$clusters, dist(X_sc), centrotypes="centroids")$DB
	coords$davies_bouldin_spectral_data  <- clusterSim::index.DB(t(metacells_integrated_data), tmp$clusters, DIST_METACELLS_INTEGRATED_DATA, centrotypes="centroids")$DB


	plotter <- reshape2::melt(coords[, c('silhouette_louvain_data', 'silhouette_Spectral_eigen', 'silhouette_Spectral_data', 'dunn_louvain_data', 'dunn_spectral_eigen','dunn_spectral_data', 'davies_bouldin_louvain_data', 'davies_bouldin_spectral_eigen', 'davies_bouldin_spectral_data')])
	plotter$method <- ifelse(grepl('dunn', plotter$variable), 'dunn', 
					ifelse(grepl('silhouette',  plotter$variable), 'silhouette', 'davies_bouldin'))
	title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Number of Nearest\nNeighbors: ", nn), fontface='bold')

	print(cowplot::plot_grid(
		cowplot::plot_grid(
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=k_11)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_manual(values=palette) + theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Spectral Clustering'),
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_Spectral_eigen)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Spectral silhouette eigen'), 
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_Spectral_data)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Spectral silhouette data'), 
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=dunn_spectral_eigen)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Spectral dunn eigen'), 	
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=dunn_spectral_data)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Spectral dunn data'), 

			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=integrated_snn_res.0.5)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_manual(values=palette)+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Louvain Clustering'), 
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=silhouette_louvain_data)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Louvain silhouette'), 
			ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=dunn_louvain_data)) + geom_point(size=0.5, alpha=0.8) + theme_void() +
			scale_color_distiller(palette = "Spectral")+ theme(legend.position="none", plot.title = element_text(size=6)) + ggtitle('Louvain dunn'), title, NULL,
		ncol=5), 

	cowplot::plot_grid(ggplot(plotter, aes(x= variable, y =value, fill=variable)) + 
	geom_boxplot() + theme_classic() + scale_fill_manual(values=palette) + theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~method, scale='free'), NULL, ncol=2, rel_widths=c(1, 0.3)),
	nrow=2, rel_heights=c(2, 1))
	)
}
dev.off()


# pdf(paste0(PLOT_PATH, '/Test.pdf'))
# palette_named <- setNames(palette, levels(combined_obj))
#  SCpubr::do_DimPlot(sample = combined_obj, label = TRUE, 
# 		colors.use =palette_named, split.by = "seurat_clusters", 
# 		ncol = 3,legend.position = "none")
# dev.off()



