
library(Seurat)
library(ggplot2)
library(enrichR)


get_enrichr_results <- function(data_set){
	all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
	for (CT in names(data_set)){
		message(CT)
		tmp <- data_set[[CT]]
		if (!'gene' %in% colnames(tmp)){
			tmp$gene <- rownames(tmp)
		}
		if(!'avg_logFC' %in% colnames(tmp)){
			tmp$avg_logFC <- tmp$avg_log2FC
		}
		tmp_pos <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >0,'gene']
		tmp_neg <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC <0,'gene']
		tmp_pos <- tmp[tmp$avg_logFC >0,'gene']
		tmp_neg <- tmp[tmp$avg_logFC <0,'gene']
		enriched <- enrichr(tmp_pos, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.5,]
		if(nrow(enriched) > 1){
				all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
		
		}
		enriched <- enrichr(tmp_neg, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.5,]
		if(nrow(enriched) > 1){

				all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))

		}
	}
	return(all_results)
}


get_dotplot <- function(data){
ggplot() + 
geom_point(data, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(name = 'Combined\nScore', low =scales::muted("blue"),mid = "white",high = scales::muted("red"),midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + 
scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(5, 1), breaks=c(0.01, 0.05, 0.1, 0.5), limits=c(0, 0.5)) +
theme_minimal() +   
theme(legend.position='bottom',
	axis.text.y=element_text(size=6),
	axis.text.x=element_text(size = 6, angle=45, vjust = 1, hjust = 1),
	legend.key.size = unit(0.2, "cm"), 
	legend.title=element_text(size=7), 
	axis.title=element_blank(), 
	legend.text=element_text(size=5)) 
}



get_volcano <- function(results, CType, sub= "A Vs B"){
	p <-  EnhancedVolcano::EnhancedVolcano(results,
				lab = results$gene_name,
				x = 'avg_logFC', captionLabSize=7,
				y = 'p_val_adj',col=c('black', 'black', 'black', 'red3'),
				subtitle =sub, axisLabSize = 5, title= CType,
				pCutoff = 0.05, labSize=2,
				FCcutoff = 0.2) +
				theme(legend.position='none', text = element_text(family = "Helvetica", size = 2),
				axis.title= element_text(family = "Helvetica", size = 7))
	return(p)
}

get_tornado <- function(results, title){
	p <- ggplot(results, aes(x = NES, y = reorder(pathway, NES), fill = is_pos)) + 
			geom_bar(stat = "identity", position = "identity") + theme_classic() +
			ylab('Pathways in C2')+
			theme(legend.position="none", text = element_text(family = "Helvetica", size = 7), axis.text.y= element_text(size = 5),) + 
			scale_fill_manual(values=c(Neg='#006e90', Pos='#cc2936')) + ggtitle(title)
	return(p)
}


scanpy_colors <- c(
		"#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
		"#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
		"#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
		"#f3e1eb", "#f6c4e1", "#f79cd4", "#7f7f7f", "#c7c7c7", "#1CE6FF", "#336600"
)

scanpy_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')


PLOT_PATH = '/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/'
combined_obj <- readRDS('/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated.rds')
scanpy_colors <- setNames(scanpy_colors, unique(combined_obj$integrated_snn_res.0.5))


# set the new clinical classification
Secondary_samples = c('aml2', 'aml7', 'aml9', 'aml11')

coords <- FetchData(combined_obj, vars = c(c('UMAP_1', 'UMAP_2', 'Sample', 'integrated_snn_res.0.5')))
coords$clinical_classification <- ifelse(coords$Sample %in% Secondary_samples, 'Secondary', 'Primary')



plot_list <- list()
for (smp in unique(coords$Sample)) {

	p <- ggplot() +
	geom_point(data = coords[coords$Sample != smp,], mapping = aes(x = UMAP_1, y = UMAP_2),  size = 0.2, color='gray', alpha=0.4) +
	geom_point(data = coords[coords$Sample == smp,], mapping = aes(x = UMAP_1, y = UMAP_2, fill = integrated_snn_res.0.5), shape=21) +
	  geom_point(size = 0.5) + scale_fill_manual(values = scanpy_colors) +
	  theme_void() + theme(legend.position = 'none') +
	  labs(color = 'Louvain clustering') +
	  ggtitle(smp)
	plot_list[[smp]] <- p
}


pdf(paste0(PLOT_PATH, 'PerSample_UMAP_and_clinical_classification.pdf'))
cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggplot(coords, aes(x = UMAP_1, y = UMAP_2, color = clinical_classification)) +
  geom_point(size = 0.5) + scale_color_manual(values = c('#AEF098', '#962493')) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(color = 'Clinical classification')

coords_tmp <- coords
coords_tmp[,'integrated_snn_res.0.5'] <- NULL
ggplot() +
geom_point(data = coords_tmp, mapping = aes(x = UMAP_1, y = UMAP_2), alpha=0.6, color='gray') +
geom_point(data = coords, mapping = aes(x = UMAP_1, y = UMAP_2, fill = integrated_snn_res.0.5), shape=21) +
	geom_point(size = 0.5) + scale_fill_manual(values = scanpy_colors) +
	theme_void() + theme(legend.position = 'none') +
	labs(color = 'Louvain clustering') + facet_wrap(~ integrated_snn_res.0.5)


coords_tmp <- coords
coords_tmp[,'integrated_snn_res.0.5'] <- NULL
ggplot() +
geom_point(data = coords_tmp, mapping = aes(x = UMAP_1, y = UMAP_2), alpha=0.6, color='gray') +
geom_point(data = coords, mapping = aes(x = UMAP_1, y = UMAP_2, fill = integrated_snn_res.0.5), shape=21) +
	geom_point(size = 0.5) + scale_fill_manual(values = scanpy_colors) +
	theme_void() + theme(legend.position = 'none') +
	labs(color = 'Louvain clustering') + facet_wrap(Sample~ integrated_snn_res.0.5)


dev.off()

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human', 'Reactome_2022')



markers <- list()
for (cl in unique(combined_obj$integrated_snn_res.0.5)){
	markers[[cl]] <- FindMarkers(combined_obj, ident.1 = cl, min.pct = 0.25, logfc.threshold = 0.25)
}

enrich_results <- get_enrichr_results(markers)
# p <- get_dotplot(enrich_results)

enrich_results[enrich_results$Combined.Score > 200, 'Combined.Score'] <- 200
enrich_results[enrich_results$Combined.Score < -200, 'Combined.Score'] <- -200

enrich_results$GO_term <- stringr::str_extract(enrich_results$Term, '(GO:[\\d]+)')

pdf(paste0(PLOT_PATH, 'ORA_per_Cluster.pdf'), width=8, height=11)
for (ct in sort(unique(enrich_results$cell_type))){
	print(ct)
	p <- get_dotplot(enrich_results[enrich_results$cell_type == ct,])
	p <- p + labs(title = ct)
	print(p)
}
dev.off()


enrich_results_list <- split(enrich_results, enrich_results$cell_type)
WriteXLS::WriteXLS(enrich_results_list, SheetNames=names(enrich_results_list), paste0(PLOT_PATH,'/Enrichr_ALL.xlsx'), row.names=TRUE)
WriteXLS::WriteXLS(markers, SheetNames=names(markers), paste0(PLOT_PATH,'/DE_ALL.xlsx'), row.names=TRUE)


pathways_c2 <- fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/c2.cp.v7.5.1.symbols.gmt')



pdf(paste0(PLOT_PATH, '/Volcano_And_tornado_1VsAll.pdf'))
for (cl in sort(names(markers))){
	print(cl)
	tmp <- markers[[cl]]
	tmp$gene <- rownames(tmp)
	tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_log2FC) )
	print(get_volcano(tt, cl, sub=paste0('Cluster ', cl, ' Vs All')))
	tmp <- tmp[order(tmp$avg_log2FC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_log2FC, tmp$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
	print(get_tornado(results_c2, cl))
	}
}
dev.off()


# Get some celltype scores and DE tests

marker_list <- list()
marker_list[['hsc']]       <-c('CRHBP', 'HOPX', 'KYT', 'CD34', 'AVP', 'PRSS2','MLLT3', 'IDS', 'BTS2')  
marker_list[['lmpp']]      <-c('PTPRC', 'FLT3', 'PROM1', 'SATB1', 'ACY3','IGJ', 'LSP1', 'NEGR1', 'SPINK2')  
marker_list[['gmp']]       <-c('CSF3R', 'CTSG', 'PRTN3', 'MPO','CFD', 'CTSG', 'CSTA', 'CST7')  
marker_list[['granul']]    <-c('ELANE', 'AZU1', 'CEBPA', 'CEBPE', 'CST7','RNASE2')  
marker_list[['mono']]      <-c('LYZ', 'CSTA', 'FCER1G', 'TLR1', 'IFNGR2')  
marker_list[['dc']]        <-c('IRF8', 'IRF7', 'IL3RA', 'CLEC4', 'IRF4', 'ZEB2', 'KLF4', 'AXL')  
marker_list[['t']]         <-c('JCHAIN', 'IKZF1', 'CYTH1', 'GZMA', 'SH2D1A', 'IL32', 'MT1E', 'CXCR3')  
marker_list[['clp']]       <-c('IL7R', 'DNTT', 'CYGB', 'LTB', 'IGJ', 'DNTT', 'ADA')  
marker_list[['prob']]      <-c('VPREB1', 'EBF1', 'CD79A', 'CD79B', 'TCL1A', 'IRF4', 'CCDC42B', 'FCRLA', 'IGLL5')  
marker_list[['mep']]       <-c('NFE2', 'HFS1', 'TAL1', 'FCER1A', 'PBX1', 'PLEK', 'DAD1', 'IGSF10')  
marker_list[['mk']]        <-c('PBX1', 'MPL', 'VWF', 'FLI1', 'ITGA22B', 'GP1BA', 'CMTM5', 'PF4', 'GP9', 'CLEC1B', 'PPBP')  
marker_list[['ery_early']] <- c('CNRIP1', 'SLC40A1', 'PKIG', 'PDZD8', 'FCER1A', 'CSF2RB', 'EMP3', 'HBD')    
marker_list[['ery_late']]  <-c('APOC1', 'BLVRB', 'APOE', 'FAM178B', 'CA1', 'CNRIP1', 'AHSP', 'KCNH2', 'TFR2', 'HBB', 'KLF1', 'S100A6')  
marker_list[['baso']]      <-c('RUNX1', 'HDC', 'MS4A2', 'MS4A3', 'TPSAB1')

rank <- 11
res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
combined_obj$nmf_clustering_11 <- clustering_nmf

DefaultAssay(combined_obj) <- 'RNA'
for (celltype in names(marker_list)){
	message(celltype)
	features <- list(c(marker_list[[celltype]]))
	tmp <- AddModuleScore(combined_obj, features=features, name = paste0(celltype, '_score'))
	coords <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.0.5', 'nmf_clustering_11',paste0(celltype, '_score1')))
	coords <- coords[order(coords[[paste0(celltype, '_score1')]]),]
	plot <- ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= get(paste0(celltype, '_score1')))) + 
			geom_point(size=1, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = celltype) + ggtitle(paste0(celltype, '_score1')) + 
			theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank())
	pdf(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/UMAP_Cluster_',celltype,'_Score_NMF.pdf'))
		print(plot)
		print(plot + 
			geom_point(data = transform(coords, nmf_clustering_11 = NULL), colour = "grey85", size=0.2, alpha=0.3) +
			geom_point(size=1, alpha=0.6) + facet_wrap(~nmf_clustering_11))
			# + facet_wrap(~integrated_snn_res.0.5))
	dev.off()
}




# Get number of cells per Seacell

library(Seurat)

seacells_origin <- list()
folders <- list.dirs('/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells' )
for (fld in folders){
	aml_name <- stringr::str_extract(fld, '(?<=\\/)aml[\\d]{1,2}$')
	if(!is.na(aml_name)){
		# read the csv
		cell_2_seacells <- read.csv(paste0(fld, '/AML_cells_2_SEACells.csv'), sep='\t', header=TRUE)
		cell_2_seacells$Sample <- aml_name
		seacells_origin[[aml_name]] <- cell_2_seacells
	}
}

seacells_origin <- do.call("rbind", seacells_origin)

seacells_origin$SEACell <- factor(seacells_origin$SEACell, levels=rev(gtools::mixedsort(unique(seacells_origin$SEACell))))
seacells_origin$Sample <- factor(seacells_origin$Sample, levels=gtools::mixedsort(unique(as.character(seacells_origin$Sample))))

seacells_origin$SEACell_name <- paste0(seacells_origin$Sample, '_', seacells_origin$SEACell)

combined_obj <- readRDS('/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated.rds')
cluster_names <- FetchData(combined_obj, vars='seurat_clusters')

seacells_origin <- merge(seacells_origin,cluster_names, by.x='SEACell_name', by.y=0)

seacells_origin_list <- split(seacells_origin, seacells_origin$Sample)

tmp <- list()
for (nm in names(seacells_origin_list)){
	tmp[[nm]] <-  as.data.frame.matrix(table(seacells_origin_list[[nm]]$seurat_clusters, seacells_origin_list[[nm]]$SEACell))
}

# write the resuslting table to an excell
WriteXLS::WriteXLS(tmp, SheetNames=c(names(tmp)), '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells_origin.xlsx', row.names=TRUE)

a <- as.data.frame(table(seacells_origin$SEACell_name))
a <- merge(as.data.frame(colnames(combined_obj)), a, by.x=1, by.y='Var1', all.x=TRUE)
rownames(a) <- a[,1]
combined_obj$cells_per_SEACell <- a[,2]

Secondary_samples = c('aml2', 'aml7', 'aml9', 'aml11')
combined_obj$clinical_classification  <- ifelse(combined_obj$Sample %in% Secondary_samples, 'Secondary', 'Primary')

saveRDS(combined_obj, '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')


coords <- FetchData(combined_obj, vars = c('UMAP_1', 'UMAP_2','cells_per_SEACell', 'Sample', 'clinical_classification'))

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/UMAP_Clinical_Classification_CellsPerSEACElls.pdf')
	ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= log10(cells_per_SEACell))) + 
	geom_point(size=1, alpha=0.6) +
	scale_color_distiller(palette = "Spectral", direction = -1, name = 'log10 Cells per SEACell') + 
	theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
	axis.ticks=element_blank(), axis.text=element_blank()) 
	ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= cells_per_SEACell)) + 
	geom_point(size=1, alpha=0.6) +
	scale_color_distiller(palette = "Spectral", direction = -1, name = 'Cells per SEACell') + 
	theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
	axis.ticks=element_blank(), axis.text=element_blank()) +facet_wrap(~Sample)
dev.off()



# DE analysis per cluster for Primary Vs Secondary
DefaultAssay(combined_obj) <- 'RNA'
results_DE <- list()
metadata <- FetchData(combined_obj, vars=c('Sample',  'clinical_classification', 'integrated_snn_res.0.5'))
combined_obj_tmp <- combined_obj
for(cl in sort(unique(combined_obj$integrated_snn_res.0.5))){
	cell_1 <-  rownames(metadata[metadata$integrated_snn_res.0.5 == cl & metadata$clinical_classification == 'Primary',])
	cell_2 <-  rownames(metadata[metadata$integrated_snn_res.0.5 == cl & metadata$clinical_classification == 'Secondary',])
	if(length(cell_1) == 0 || length(cell_2) == 0){
		next
	}
	combined_obj_tmp$Status <- ifelse(colnames(combined_obj_tmp) %in% cell_1, 'Primary', 
									 ifelse(colnames(combined_obj_tmp) %in% cell_2, 'Secondary', 'None'))
	Idents(combined_obj_tmp)  <- 'Status'
	results <- FindMarkers(combined_obj_tmp, 
								ident.1 = 'Primary',
								ident.2 = 'Secondary',
								assay='RNA',
								test.use = 'MAST', 
								latent.vars='Sample')
	results_DE[[cl]] <- results
}
WriteXLS::WriteXLS(results_DE, SheetNames=names(results_DE), '/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/DE_per_Cluster_Primary_Vs_Secondary.xlsx', row.names=TRUE)



DefaultAssay(combined_obj) <- 'RNA'
Idents(combined_obj) <- 'integrated_snn_res.0.5'
list_DE <- list()
for(cl1 in sort(unique(combined_obj$integrated_snn_res.0.5))){
	for (cl2 in sort(unique(combined_obj$integrated_snn_res.0.5))){
		if (as.numeric(cl1) >= as.numeric(cl2)){
			next
		}
		else{
			print(paste0(cl1, ' vs ', cl2))
			list_DE[[paste0(cl1, 'vs', cl2)]]<- FindMarkers(
			object = combined_obj,
			ident.1 = cl1,
			ident.2 = cl2,
			assay = "RNA",
			recorrect_umi = FALSE)
		}
	}
}
WriteXLS::WriteXLS(list_DE, ExcelFileName='/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/DE_Clusters_1Vs_1.xlsx', SheetNames = names(list_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)






combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')

# Plot the PCA for each SEACell and color it by sample or condition
coords <- FetchData(combined_obj, vars = c('PC_1', 'PC_2','Sample', 'clinical_classification', 'integrated_snn_res.0.5'))

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/PCA_sample_clinic.pdf')
	ggplot(coords, aes(x=PC_1, y= PC_2, color= Sample)) + 
		geom_point(size=1, alpha=0.8) +
		scale_color_manual(values=scanpy_colors) + 
		theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
		axis.ticks=element_blank(), axis.text=element_blank())
	ggplot(coords, aes(x=PC_1, y= PC_2, color= clinical_classification)) + 
		geom_point(size=1, alpha=0.8) +
		scale_color_manual(values=scanpy_colors) +
		theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
		axis.ticks=element_blank(), axis.text=element_blank())
dev.off()


# Cells per metacell

plotter <- FetchData(combined_obj, vars = c('cells_per_SEACell', 'Sample', 'clinical_classification'))
plotter$SEACell <- stringr::str_extract(rownames(plotter), 'SEACell-[\\d]+')

plotter$SEACell <- factor(plotter$SEACell, levels=gtools::mixedsort(unique(plotter$SEACell), decreasing=TRUE))

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/Cells_per_Metacell_perSample.pdf')
	ggplot(plotter, aes(x=SEACell, y=cells_per_SEACell, group=Sample, color= Sample)) + 
		geom_line(alpha=0.8) +
		scale_color_manual(values=scanpy_colors) + 
		theme_classic() + theme(legend.position="bottom", panel.background = element_blank(), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size = 6)) + 
		guides(color=guide_legend(nrow=2,byrow=TRUE))
		ggplot(plotter, aes(x=SEACell, y=cells_per_SEACell, group=Sample, color= Sample)) + 
		geom_line(alpha=0.8) +
		scale_color_manual(values=scanpy_colors) + 
		theme_classic() + theme(legend.position="none", panel.background = element_blank(), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
		guides(color=guide_legend(nrow=2,byrow=TRUE)) + facet_wrap(~Sample, ncol=3)
dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/Cells_per_Metacell_clinical_classification.pdf')
	ggplot(plotter, aes(x=SEACell, y=cells_per_SEACell, group=clinical_classification, color= clinical_classification)) + 
		geom_line(alpha=0.8) +
		scale_color_manual(values=scanpy_colors) + 
		theme_classic() + theme(legend.position="bottom", panel.background = element_blank(), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size = 6)) + 
		guides(color=guide_legend(nrow=2,byrow=TRUE))
		ggplot(plotter, aes(x=SEACell, y=cells_per_SEACell, group=clinical_classification, color= clinical_classification)) + 
		geom_line(alpha=0.8) +
		scale_color_manual(values=scanpy_colors) + 
		theme_classic() + theme(legend.position="none", panel.background = element_blank(), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
		guides(color=guide_legend(nrow=2,byrow=TRUE)) + facet_wrap(~clinical_classification, nrow=2)
dev.off()




# Signatures from C2 and C5 msigdb



C5_C2_pathways <- c(
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/c5.all.v2023.1.Hs.symbols.gmt'), 
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/c2.all.v2023.1.Hs.symbols.gmt'), 
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/GOBP_RESPONSE_TO_DRUG.v2023.1.Hs.gmt')
	)

C5_terms <- c(
	'GOBP_HEMATOPOIETIC_STEM_CELL_PROLIFERATION',
	'GOBP_STEM_CELL_PROLIFERATION',
	'GOBP_POSITIVE_REGULATION_OF_HEMATOPOIETIC_STEM_CELL_PROLIFERATION',
	# 'GOBP_NEGATIVE_REGULATION_OF_HEMATOPOIETIC_STEM_CELL_PROLIFERATION',
	'GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION',
	'GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION', #'G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE',
	'GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION', #'G2_M_TRANSITION_OF_MITOTIC_CELL_CYCLE',
	'GOBP_CELL_CYCLE_DNA_REPLICATION',
	'GOBP_NEGATIVE_REGULATION_OF_CELL_CYCLE',
	'GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PROCESS',
	'KEGG_CELL_CYCLE'
	# 'M_PHASE_OF_MITOTIC_CELL_CYCLE'
)

C2_terms <- c(
	'KESHELAVA_MULTIPLE_DRUG_RESISTANCE',
	'GOBP_REGULATION_OF_RESPONSE_TO_DRUG',
	'GOBP_RESPONSE_TO_DRUG',
	'YAGUE_PRETUMOR_DRUG_RESISTANCE_DN',
	'YAGUE_PRETUMOR_DRUG_RESISTANCE_UP',
	'WP_LNCRNAMEDIATED_MECHANISMS_OF_THERAPEUTIC_RESISTANCE',
	'TOOKER_GEMCITABINE_RESISTANCE_DN',
	'TOOKER_GEMCITABINE_RESISTANCE_UP',
	'TSUNODA_CISPLATIN_RESISTANCE_DN',
	'TSUNODA_CISPLATIN_RESISTANCE_UP',
	'WHITESIDE_CISPLATIN_RESISTANCE_DN',
	'WHITESIDE_CISPLATIN_RESISTANCE_UP',
	'WP_EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE',
	'RIGGINS_TAMOXIFEN_RESISTANCE_DN',
	'RIGGINS_TAMOXIFEN_RESISTANCE_UP',
	'MASSARWEH_TAMOXIFEN_RESISTANCE_DN',
	'MASSARWEH_TAMOXIFEN_RESISTANCE_UP',
	'MCCOLLUM_GELDANAMYCIN_RESISTANCE_DN',
	'MCCOLLUM_GELDANAMYCIN_RESISTANCE_UP',
	'LI_CISPLATIN_RESISTANCE_DN',
	'LI_CISPLATIN_RESISTANCE_UP',
	'MAHADEVAN_IMATINIB_RESISTANCE_DN',
	'MAHADEVAN_IMATINIB_RESISTANCE_UP',
	'MARCHINI_TRABECTEDIN_RESISTANCE_DN',
	'MARCHINI_TRABECTEDIN_RESISTANCE_UP',
	'MASRI_RESISTANCE_TO_TAMOXIFEN_AND_AROMATASE_INHIBITORS_DN',
	'MASRI_RESISTANCE_TO_TAMOXIFEN_AND_AROMATASE_INHIBITORS_UP',
	'HUANG_DASATINIB_RESISTANCE_DN',
	'KANG_CISPLATIN_RESISTANCE_DN',
	'KANG_CISPLATIN_RESISTANCE_UP',
	'KANG_DOXORUBICIN_RESISTANCE_DN',
	'KANG_DOXORUBICIN_RESISTANCE_UP',
	'KANG_FLUOROURACIL_RESISTANCE_DN',
	'KANG_FLUOROURACIL_RESISTANCE_UP',
	'HOLLEMAN_VINCRISTINE_RESISTANCE_ALL_DN',
	'HOLLEMAN_VINCRISTINE_RESISTANCE_ALL_UP',
	'HOLLEMAN_PREDNISOLONE_RESISTANCE_ALL_DN',
	'HOLLEMAN_PREDNISOLONE_RESISTANCE_ALL_UP',
	'GYORFFY_DOXORUBICIN_RESISTANCE',
	'GYORFFY_MITOXANTRONE_RESISTANCE',
	'HANN_RESISTANCE_TO_BCL2_INHIBITOR_DN',
	'HANN_RESISTANCE_TO_BCL2_INHIBITOR_UP',
	'HOLLEMAN_ASPARAGINASE_RESISTANCE_ALL_DN',
	'HOLLEMAN_ASPARAGINASE_RESISTANCE_ALL_UP'
)

DefaultAssay(combined_obj) <- 'RNA'
combined_obj <- ScaleData(combined_obj)
all_terms <- c(C5_terms, C2_terms)

all_signatures <- list()
for (term in all_terms){
	print(term)
	features <- list(C5_C2_pathways[[term]])
	tmp <- FetchData(AddModuleScore(combined_obj, features=features, name = term), 
					vars = c('UMAP_1', 'UMAP_2', 'Sample', paste0(term, '1')))
	tmp$cell_id <- rownames(tmp)
	all_signatures[[term]] <- tmp
	# if(nrow(all_signatures) ==0){
	# 	all_signatures <- tmp
	# }else{
	# 	all_signatures <- merge(all_signatures, tmp[, c('cell_id', paste0(term, '1'))], by='cell_id')
	# }
}

plot_signature_UMAP <- function(signature_data){
	signature_data <- signature_data[order(signature_data[[4]], decreasing=FALSE),]
	p <- ggplot(signature_data, aes(x=UMAP_1, y=UMAP_2, color=signature_data[[4]])) + 
		geom_point(size=1, alpha=0.8) +
		ggtitle(sub('1$', '', names(signature_data)[4], perl=TRUE)) +
		scale_color_distiller(palette = "Spectral", direction = -1) + 
		theme_void() + theme(
			legend.position="bottom", 
			panel.background = element_blank(), 
			axis.ticks=element_blank(), 
			axis.text=element_blank(),
			legend.title=element_blank(), 
			plot.title=element_text(size=6)) +
		guides(colour=guide_colourbar(barwidth=10, barheight=0.5))

	return(p)
}


pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/Signatures.pdf')
for(i in seq(from = 1, to = length(all_signatures), by = 4)){
	plotter <- list()
	plotter[[1]] <- plot_signature_UMAP(all_signatures[[i]])
	if(i+1 < length(all_signatures)){
		plotter[[2]] <- plot_signature_UMAP(all_signatures[[i+1]])
	}
	if(i+2 < length(all_signatures)){
		plotter[[3]] <- plot_signature_UMAP(all_signatures[[i+2]])
	}
	if(i+3 < length(all_signatures)){
		plotter[[4]] <- plot_signature_UMAP(all_signatures[[i+3]])
	}
	print(cowplot::plot_grid(plotlist = plotter, ncol=2, nrow=2))
}
dev.off()


# NMF clustering


library(Seurat)
library(ggplot2)

combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')

# data_counts <- as.matrix(combined_obj@assays$RNA@counts)
rank.range <- seq(5, 20)
# for (rank in rank.range){
# 	print(rank)
# 	res <- NMF::nmf(data_counts, rank)
# 	saveRDS(res, paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_',rank,'.rds'))
# }



clustering_results <- list()
for (rank in rank.range){
	print(rank)
	# res <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_',rank,'.rds'))
	# res <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Norm_',rank,'.rds'))
	res <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
	clustering_results[[rank]] <- res
}

# get the silhouette.coef from the NMF clustering
plotter <- data.frame(Rank=numeric(), Silhouette=numeric())
for (rank in rank.range){
	print(rank)
	res <- clustering_results[[rank]]
	plotter <- rbind(plotter, data.frame(Rank=rank, Silhouette=NMF::summary(res)[['silhouette.coef']]))
}


scanpy_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', 
'#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')


coords <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'Sample'))

for (rank in rank.range){
	coords <- merge(coords, setNames(as.data.frame(NMF::predict(clustering_results[[rank]])), paste0('Cluster_', rank)), by=0)
	rownames(coords) <- coords[,'Row.names']
	coords <- coords[, !names(coords) %in% 'Row.names']
}


get_umap_rank <- function(data, rank){
	cluster_rank <- paste0('Cluster_', rank)
	if(length(unique(data[[cluster_rank]]))>16){
		scanpy_colors_tmp <- colorRampPalette(scanpy_colors)(length(unique(data[[cluster_rank]])))
	}else{
		scanpy_colors_tmp <- scanpy_colors
	}
	p <- ggplot() +
		geom_point(data = data, mapping = aes(x = UMAP_1, y = UMAP_2, fill = get(cluster_rank)), shape=21, size = 1.5) +
		scale_fill_manual(values = scanpy_colors_tmp) +
		theme_void() + theme(legend.position = 'none') +
		ggtitle(paste0('NMF rank: ', rank))
	return(p)
}

get_umap_faceted <-  function(data, rank){
	cluster_rank <- paste0('Cluster_', rank)
	if(length(unique(data[[cluster_rank]]))>16){
		scanpy_colors_tmp <- colorRampPalette(scanpy_colors)(length(unique(data[[cluster_rank]])))
	}else{
		scanpy_colors_tmp <- scanpy_colors
	}
	p <- ggplot() +
		geom_point(data = coords[,  c('UMAP_1', 'UMAP_2')], mapping = aes(x = UMAP_1, y = UMAP_2), alpha=0.6, color='gray', size = 0.3) +
		geom_point(data = data, mapping = aes(x = UMAP_1, y = UMAP_2, fill = get(cluster_rank)), shape=21, size = 1) +
		scale_fill_manual(values = scanpy_colors_tmp) +
		theme_void() + theme(legend.position = 'none') +
		facet_wrap(~get(cluster_rank))
	return(p)
}

# pdf('./Plots/NMF_checks/NMF_Silhouette_RawCounts.pdf')
pdf('./Plots/NMF_checks/NMF_Silhouette_Norm_Integrated.pdf')
ggplot(plotter, aes(x=Rank, y=Silhouette)) + geom_point() + geom_line() + theme_classic()
cowplot::plot_grid(get_umap_rank(coords, 5),  get_umap_faceted(coords, 5))
cowplot::plot_grid(get_umap_rank(coords, 6),  get_umap_faceted(coords, 6))
cowplot::plot_grid(get_umap_rank(coords, 7),  get_umap_faceted(coords, 7))
cowplot::plot_grid(get_umap_rank(coords, 8),  get_umap_faceted(coords, 8))
cowplot::plot_grid(get_umap_rank(coords, 9),  get_umap_faceted(coords, 9))
cowplot::plot_grid(get_umap_rank(coords, 10), get_umap_faceted(coords, 10))
cowplot::plot_grid(get_umap_rank(coords, 11), get_umap_faceted(coords, 11))
cowplot::plot_grid(get_umap_rank(coords, 12), get_umap_faceted(coords, 12))
cowplot::plot_grid(get_umap_rank(coords, 13), get_umap_faceted(coords, 13))
cowplot::plot_grid(get_umap_rank(coords, 14), get_umap_faceted(coords, 14))
cowplot::plot_grid(get_umap_rank(coords, 15), get_umap_faceted(coords, 15))
cowplot::plot_grid(get_umap_rank(coords, 16), get_umap_faceted(coords, 16))
cowplot::plot_grid(get_umap_rank(coords, 17), get_umap_faceted(coords, 17))
cowplot::plot_grid(get_umap_rank(coords, 18), get_umap_faceted(coords, 18))
cowplot::plot_grid(get_umap_rank(coords, 19), get_umap_faceted(coords, 19))
cowplot::plot_grid(get_umap_rank(coords, 20), get_umap_faceted(coords, 20))
dev.off()

pdf('./Plots/NMF_checks/NMF_Silhouette_Norm_PerSample_Integrated.pdf')
get_umap_rank(coords, 20)
ggplot(coords, aes(x = UMAP_1, y = UMAP_2, fill = Cluster_20)) + 
	geom_point(shape=21, size = 1.5) +
	scale_fill_manual(values = colorRampPalette(scanpy_colors)(20)) +
	theme_void() + theme(legend.position = 'none') +
	facet_wrap(~Sample)
ggplot() + 
	geom_point(data = coords[, c('UMAP_1','UMAP_2', 'Sample')], mapping=aes(x = UMAP_1, y = UMAP_2), color='gray', alpha=0.6,  size = 0.5) +
	geom_point(data = coords, mapping=aes(x = UMAP_1, y = UMAP_2, fill = Cluster_20), shape=21, size = 1.5) +
	scale_fill_manual(values = colorRampPalette(scanpy_colors)(20)) +
	theme_void() + theme(legend.position = 'none') +
	facet_wrap(~Cluster_20)
dev.off()


res <- clustering_results[[15]]
pdf('./Plots/NMF_checks/NMF_Silhouette_coefs.pdf')
NMF::consensusmap(res)
NMF::basismap(res, subsetRow=TRUE)
NMF::coefmap(res)
dev.off()





# NMF rank: 17

combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')
rank <- 7
res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
combined_obj$nmf_clustering_7 <- clustering_nmf
rank <- 11
res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
combined_obj$nmf_clustering_11 <- clustering_nmf
rank <- 14
res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
combined_obj$nmf_clustering_14 <- clustering_nmf

coords <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'Sample', 'nmf_clustering_7', 'nmf_clustering_11', 'nmf_clustering_14', 'integrated_snn_res.0.5'))


scanpy_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', 
'#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')
scanpy_colors_tmp <- colorRampPalette(scanpy_colors)(17)


get_umap <- function(data, color_by){
	if(color_by == 'integrated_snn_res.0.5'){
		title <- 'Louvain'
	} else if(stringr::str_detect(color_by, 'nmf')){
		title <- paste0('NMF ', stringr::str_extract(color_by, '[0-9]+'))
	} else{
		title <- color_by
	}
	p <- ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=coords[[color_by]])) + 
		 geom_point(size=0.5, alpha=1) + scale_color_manual(values=scanpy_colors_tmp) +
		 theme_void() + theme(legend.position="none") + ggtitle(title)
	return(p)
}


coords$Sample <- factor(coords$Sample, levels =gtools::mixedsort(unique(coords$Sample)))

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/RPCA_integrated/NMF_Vs_Louvain.pdf')
cowplot::plot_grid(
	cowplot::plot_grid(
		NULL, 
		get_umap(coords, 'integrated_snn_res.0.5'),
		NULL, 
	ncol=3, rel_widths=c(0.25,0.5,0.25))
	,	
	cowplot::plot_grid(
		get_umap(coords, 'nmf_clustering_7'),
		get_umap(coords, 'nmf_clustering_11'),
		get_umap(coords, 'nmf_clustering_14'),
	ncol=3)
,nrow=2)

ggplot() + 
	geom_point(data = coords[, c('UMAP_1','UMAP_2')], mapping=aes(x = UMAP_1, y = UMAP_2), shape=16, color='gray', alpha=0.6,  size = 0.5) +
	geom_point(data = coords, mapping=aes(x = UMAP_1, y = UMAP_2, color = nmf_clustering_11), size = 1.5) +
	scale_color_manual(values = scanpy_colors_tmp) +
	theme_void() + theme(legend.position = 'none') +
	facet_wrap(~Sample) + ggtitle('NMF 11 per sample')
dev.off()


# We are keeping the NMF 11