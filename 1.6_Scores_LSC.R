

library(Seurat)
library(ggplot2)

scanpy_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', 
'#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')
scanpy_colors_tmp <- colorRampPalette(scanpy_colors)(17)


combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')


coords <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'Sample', 'nmf_clustering_7', 'nmf_clustering_11', 'nmf_clustering_14', 'integrated_snn_res.0.5'))

get_umap <- function(data, color_by){
	if(color_by == 'integrated_snn_res.0.5'){
		title <- 'Louvain'
	} else if(stringr::str_detect(color_by, 'nmf')){
		title <- paste0('NMF ', stringr::str_extract(color_by, '[0-9]+'))
	} else{
		title <- color_by
	}
	p <- ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=.data[[color_by]])) + 
		 geom_point(size=0.5, alpha=1) + scale_color_manual(values=scanpy_colors_tmp) +
		 theme_void() + theme(legend.position="right") + ggtitle(title)
	return(p)
}

get_umap_signature <- function(data, color_by_pos, title){
	p <- ggplot(data, aes(x=UMAP_1, y= UMAP_2, color= .data[[color_by_pos]])) + 
			geom_point(size=1, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = title) + ggtitle(title) + 
			theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank()) 
	return(p)
}

get_violin <- function(data, color_by_pos, title, clustering_res = 'nmf_clustering_11'){
	p <- ggplot(data, aes(x=.data[[clustering_res]], y=.data[[color_by_pos]], fill=nmf_clustering_11)) + 
			geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill="white")+
			scale_fill_manual(values=scanpy_colors_tmp) + labs(y=title) + 
			theme_classic() + theme(legend.position="none", panel.background = element_blank())
	return(p)
}



# Extract the MetaFeatures from the NMF clustering
metagenes_rank <- list()
for (rank in c(7, 11, 14)){
	print(rank)
	res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
	clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
	genes <- rownames(res_clustering@fit@W)
	metagenes <- NMF::extractFeatures(res_clustering, 30)
	metagenes_per_cluster <- data.frame(gene_name=NULL, cluster=NULL)
	for (clt in unique(clustering_nmf[[1]])){
		names <- genes[metagenes[[as.numeric(clt)]]]
		metagenes_per_cluster <- rbind(metagenes_per_cluster, data.frame(gene_name=names, cluster=clt))
	}
	metagenes_rank[[as.character(rank)]] <- metagenes_per_cluster
}

WriteXLS::WriteXLS(metagenes_rank, paste0('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/NMF_checks/genes_NMF_Rank_',rank,'.xlsx'), SheetNames=names(metagenes_rank))


pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/NMF_checks/genes_NMF_Rank_UMAPS.pdf')

get_umap(coords, 'nmf_clustering_7')
get_umap(coords, 'nmf_clustering_11')
get_umap(coords, 'nmf_clustering_14')

dev.off()






################################################
# CALCULATING THE DIFFERENT SCORES

# those genes not in the dataset:
# "KIAA0125" "NGFRAP1"  "GPR56" 
# https://www.nature.com/articles/nature20598
genes_4_score <- c('DNMT3B', 'ZBTB46', 'NYNRIN', 'ARHGAP22', 'LAPTM4B', 'MMRN1', 'DPYSL3',  'CDK6', 'CPXM1', 'SOCS2', 'SMIM24', 'EMP1',  'CD34', 'AKR1C3')
genes_4_score_list <- list(genes_4_score)
names(genes_4_score_list) <- 'all_genes_signature'

positive_genes <- list(c('DNMT3B',
					'NYNRIN',
					'LAPTM4B',
					'MMRN1',
					'DPYSL3',
					'SOCS2',
					'EMP1',
					'CD34'))
names(positive_genes) <- 'positive_genes_signature'
negative_genes <- list(c('ZBTB46', 
					'ARHGAP22', 
					'CDK6', 
					'CPXM1', 
					'SMIM24', 
					'AKR1C3'))
names(negative_genes) <- 'negative_genes_signature'

# https://www.nature.com/articles/nm.2415
new_gene_list <- list(c('NPAL2', 'RBPMS', 'TRAF3IP2', 'PPP1R10', 'ATP1B1', 'NF1', 'RBPMS', 'FLJ13197', 'RBPMS', 'ABCG1', 'CLN5', 'LRRC8B', 'FRMD4B', 'ZFP30', 'C17orf86', 'C16orf5', 'TGIF2', 'RABGAP1', 'PPIG', 'GPR56', 'EIF2S3', 'NAB1', 'LRRC61', 'ATP1B1', 'ZNF500', 'CSDE1', 'C2CD2', 'PAQR6', 'FAM119B', 'ARPP-19', 'SETDB1', 'ZBTB39', 'RBPMS', 'SLC9A7', 'MAP3K7', 'ARL3', 'ZNF304', 'LOC552889', 'VGLL4', 'UBR5', 'PTCD2', 'CRKRS', 'IQGAP2', 'PLCH1', 'ARFGEF1', 'MAP3K7', 'PNPLA4'))
names(new_gene_list) <- 'new_gene_list_signature'

# done in docker container
# tmp_data <- combined_obj@assays$RNA@data
# geneSets_all      <- GSEABase::GeneSet(genes_4_score, setName = 'all_genes_signature')
# geneSets_positive <- GSEABase::GeneSet(unlist(positive_genes), setName = 'positive_genes_signature')
# geneSets_negative <- GSEABase::GeneSet(unlist(negative_genes), setName = 'negative_genes_signature')
# cells_AUC <- AUCell::AUCell_run(tmp_data, geneSets_all)
# saveRDS(cells_AUC, '/SEACells_AML/Data/Other/Signature_geneSets_all.rds')
# cells_AUC <- AUCell::AUCell_run(tmp_data, geneSets_positive)
# saveRDS(cells_AUC, '/SEACells_AML/Data/Other/Signature_geneSets_posl.rds')
# cells_AUC <- AUCell::AUCell_run(tmp_data, geneSets_negative)
# saveRDS(cells_AUC, '/SEACells_AML/Data/Other/Signature_geneSets_neg.rds')
cluster_umap <- ggplot(FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11')), aes(x=UMAP_1, y= UMAP_2, color= nmf_clustering_11)) + 
	geom_point(size=1, alpha=0.9) + scale_color_manual(values=scanpy_colors_tmp)+
	guides(colour = guide_legend(override.aes = list(size=2)))+
	theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
	axis.ticks=element_blank(), axis.text=element_blank()) 

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/Signature_Nat_17_Paper.pdf')
print(cluster_umap)
LSC_score <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,genes_4_score), slot='scale.data')
LSC_score$PaperScore <- (LSC_score[,'DNMT3B'] * 0.0874) + 
						(LSC_score[,'ZBTB46'] * -0.0347) + 
						(LSC_score[,'NYNRIN'] * 0.00865) + 
						(LSC_score[,'ARHGAP22'] * -0.0138) + 
						(LSC_score[,'LAPTM4B'] * 0.00582) + 
						(LSC_score[,'MMRN1'] * 0.0258) + 
						(LSC_score[,'DPYSL3'] * 0.0284) + 
						(LSC_score[,'CDK6'] * -0.0704) + 
						(LSC_score[,'CPXM1'] * -0.0258) + 
						(LSC_score[,'SOCS2'] * 0.0271) + 
						(LSC_score[,'SMIM24'] * -0.0226) + 
						(LSC_score[,'EMP1'] * 0.0146) + 
						(LSC_score[,'CD34'] * 0.0338) + 
						(LSC_score[,'AKR1C3'] * -0.0402)
cowplot::plot_grid(get_umap_signature(LSC_score, 'PaperScore', 'PaperScore\nscaled_data'),
get_violin(LSC_score, 'PaperScore', 'PaperScore'), nrow=2)

tmp <- LSC_score
tmp$binarized_score <- ifelse(tmp$PaperScore > median(tmp$PaperScore), 'LSC', 'Other')
cowplot::plot_grid(ggplot(tmp, aes(x=UMAP_1, y=UMAP_2, color=binarized_score)) + geom_point(size=1, alpha=0.9) + 
	scale_color_manual(values=c('LSC'='#279e68', 'Other'='#aa40fc'))+
	guides(colour = guide_legend(override.aes = list(size=2)))+
	theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
	axis.ticks=element_blank(), axis.text=element_blank()) + ggtitle('Median Binarized PaperScore') ,
ggplot(tmp, aes(x=nmf_clustering_11, fill=binarized_score)) + geom_bar(position='fill') + 
	scale_fill_manual(values=c('LSC'='#279e68', 'Other'='#aa40fc'))+
	theme_classic() + theme(legend.position="right", panel.background = element_blank()) + ggtitle('Median Binarized PaperScore Distribution'),
nrow=2)
dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/Signature_Nat_17_Seurat.pdf')
print(cluster_umap)
tmp <- FetchData(
				AddModuleScore(combined_obj, features=genes_4_score_list, name = 'LSC_SeuratScore_allGenes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' , 'LSC_SeuratScore_allGenes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'LSC_SeuratScore_allGenes1', 'LSC_SeuratScore\nallGenes'),
get_violin(tmp, 'LSC_SeuratScore_allGenes1', 'LSC_SeuratScore\nallGenes'), nrow=2)

tmp <- FetchData(
				AddModuleScore(combined_obj, features=positive_genes, name = 'LSC_SeuratScore_PositiveGenes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'LSC_SeuratScore_PositiveGenes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'LSC_SeuratScore_PositiveGenes1','LSC_SeuratScore\nPositiveGenes'),
get_violin(tmp, 'LSC_SeuratScore_PositiveGenes1', 'LSC_SeuratScore\nPositiveGenes'), nrow=2)

tmp <- FetchData(
				AddModuleScore(combined_obj, features=negative_genes, name = 'LSC_SeuratScore_NegativeGenes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'LSC_SeuratScore_NegativeGenes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'LSC_SeuratScore_NegativeGenes1', 'LSC_SeuratScore\nNegativeGenes'),
get_violin(tmp, 'LSC_SeuratScore_NegativeGenes1', 'LSC_SeuratScore\nNegativeGenes'), nrow=2)
dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/Signature_Nat_17_Ucell.pdf')
print(cluster_umap)
tmp <- FetchData(
				UCell::AddModuleScore_UCell(combined_obj, features=genes_4_score_list, name='_UCell'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'all_genes_signature_UCell')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'all_genes_signature_UCell', 'LSC_UCellScore\nallGenes'),
get_violin(tmp, 'all_genes_signature_UCell', 'all_genes\nsignature_UCell'), nrow=2)

tmp <- FetchData(
				UCell::AddModuleScore_UCell(combined_obj, features=positive_genes, name='_UCell'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'positive_genes_signature_UCell')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'positive_genes_signature_UCell', 'LSC_UCellScore\nPositiveGenes'),
get_violin(tmp, 'positive_genes_signature_UCell', 'positive_genes\nsignature_UCell'), nrow=2)

tmp <- FetchData(
				UCell::AddModuleScore_UCell(combined_obj, features=negative_genes, name='_UCell'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'negative_genes_signature_UCell')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'negative_genes_signature_UCell', 'LSC_UCellScoreNegativeGenes'),
get_violin(tmp, 'negative_genes_signature_UCell', 'negative_genes\nsignature_UCell'), nrow=2)
dev.off()


pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/Signature_Nat_17_AUCell.pdf')
print(cluster_umap)
coords <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11'))
tmp <- merge(coords, 
	t(as.data.frame(AUCell::getAUC(readRDS('./Data/Other/Signature_geneSets_all.rds')))), 
	by=0)
cowplot::plot_grid(get_umap_signature(tmp, 'all_genes_signature', 'LSC_AUCellS\nallGenes'),
get_violin(tmp, 'all_genes_signature', 'all_genes_signature'), nrow=2)

tmp <- merge(coords, 
	t(as.data.frame(AUCell::getAUC(readRDS('./Data/Other/Signature_geneSets_posl.rds')))), 
	by=0)
cowplot::plot_grid(get_umap_signature(tmp, 'positive_genes_signature', 'LSC_AUCellS\npositiveGenes'),
get_violin(tmp, 'positive_genes_signature', 'LSC_AUCellSpositiveGenes'), nrow=2)

tmp <- merge(coords, 
	t(as.data.frame(AUCell::getAUC(readRDS('./Data/Other/Signature_geneSets_neg.rds')))), 
	by=0)
cowplot::plot_grid(get_umap_signature(tmp, 'negative_genes_signature', 'LSC_AUCellS\nnegativeGenes'),
get_violin(tmp, 'negative_genes_signature', 'LSC_AUCellSnegativeGenes'), nrow=2)
dev.off()


pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/NewNereaNat_Signature_Nat_UcellSeurat.pdf')
print(cluster_umap)
# new signature
tmp <- FetchData(
				AddModuleScore(combined_obj, features=new_gene_list, name = 'new_gene_list_signature_Seurat', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'new_gene_list_signature_Seurat1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'new_gene_list_signature_Seurat1', 'LSC_newSig\nScoreAllGenes'),
get_violin(tmp, 'new_gene_list_signature_Seurat1', 'signature_Seurat'), nrow=2)

tmp <- FetchData(
				UCell::AddModuleScore_UCell(combined_obj, features=new_gene_list, name='_UCell'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'new_gene_list_signature_UCell')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'new_gene_list_signature_UCell', 'LSC_UCellScoreAllGenes'),
get_violin(tmp, 'new_gene_list_signature_UCell', 'signature_UCell'), nrow=2)
dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/NewAintzBulk_Signature_Seurat.pdf')
print(cluster_umap)
# signatures from Aintzane LSC Vs Blasts
pos_aintz_genes <- read.table('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/upreg_DEGs_LSCsVsBlasts_padj0.05.txt', sep = '\t', header = TRUE)
neg_aintz_genes <- read.table('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/downreg_DEGs_LSCsVsBlasts_padj0.05.txt', sep = '\t', header = TRUE)

all_aintz_genes <- list(c(neg_aintz_genes[,1], pos_aintz_genes[, 1]))
names(all_aintz_genes) <- 'all_aintz_genes'
tmp <- FetchData(
				AddModuleScore(combined_obj, features=all_aintz_genes, name = 'all_aintz_genes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'all_aintz_genes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'all_aintz_genes1', 'LSC_AtzSig\nScoreAllGenes'),
get_violin(tmp, 'all_aintz_genes1', 'signature_Seurat'), nrow=2)

pos_aintz_genes <- list(pos_aintz_genes[,1])
names(pos_aintz_genes) <- 'pos_aintz_genes'
tmp <- FetchData(
				AddModuleScore(combined_obj, features=pos_aintz_genes, name = 'pos_aintz_genes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'pos_aintz_genes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'pos_aintz_genes1', 'LSC_AtzSigPos\nScoreAllGenes'),
get_violin(tmp, 'pos_aintz_genes1', 'signature_Seurat'), nrow=2)

neg_aintz_genes <- list(neg_aintz_genes[,1])
names(neg_aintz_genes) <- 'neg_aintz_genes'
tmp <- FetchData(
				AddModuleScore(combined_obj, features=neg_aintz_genes, name = 'neg_aintz_genes', assay='RNA'),
				vars=c('UMAP_1', 'UMAP_2', 'nmf_clustering_11' ,'neg_aintz_genes1')
				)
cowplot::plot_grid(get_umap_signature(tmp, 'neg_aintz_genes1', 'LSC_AtzNeg\nScoreAllGenes'),
get_violin(tmp, 'neg_aintz_genes1', 'signature_Seurat'), nrow=2)
ggplot(tmp, aes(x=score)) + geom_histogram(bins= 100) + theme_classic()

dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/NewAintzBulk_Signature_CART.pdf')
print(cluster_umap)
pos_aintz_genes <- setNames(read.table('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/upreg_DEGs_LSCsVsBlasts_padj0.05.txt', sep = '\t', header = TRUE), c('gene_id', 'logFC'))
neg_aintz_genes <- setNames(read.table('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/downreg_DEGs_LSCsVsBlasts_padj0.05.txt', sep = '\t', header = TRUE), c('gene_id', 'logFC'))
all_atz_genes <- rbind(pos_aintz_genes, neg_aintz_genes)
# get the genes for the signature from the dataset
tmp <- t(FetchData(combined_obj, vars = all_atz_genes$gene_id))
# remove the genes not present in the dataset
all_atz_genes <- all_atz_genes[all_atz_genes$gene_id %in% rownames(tmp), ]

logplusone <- function(x) {log(x + 0.5)}
tmp <- apply(tmp, 2, logplusone)
tmp <- scale(tmp)

tmp <- tmp * all_atz_genes$logFC
score <- base::colSums(tmp)
score <- as.data.frame(score)
tmp <- FetchData(combined_obj, vars = c('UMAP_1', 'UMAP_2', 'nmf_clustering_11', 'nmf_clustering_14', 'nmf_clustering_7'))
tmp <- merge(tmp, score, by=0)

cowplot::plot_grid(get_umap_signature(tmp, 'score', 'LSC_Atz\nScoreCARThl')+ geom_text(aes(label=ifelse(score>=0,as.character(Row.names),'')),hjust=0,vjust=0),
get_violin(tmp, 'score', 'signature_CART'), nrow=2)
ggplot(tmp) + geom_histogram(aes(x=score), bins=100) + theme_classic()

# cowplot::plot_grid(
# 	get_umap_signature(tmp, 'score', 'quantile99\nScoreCARThl')+ geom_text(aes(label=ifelse(score>=quantile(tmp$score, probs=.99),as.character(Row.names),'')),hjust=0,vjust=0),
# 	get_umap_signature(tmp, 'score', 'quantile95\nScoreCARThl')+ geom_text(aes(label=ifelse(score>=quantile(tmp$score, probs=.95),as.character(Row.names),'')),hjust=0,vjust=0), 
# 	nrow=2
# )

ggplot() +
geom_point(data = tmp[!tmp$score >= quantile(tmp$score, probs=.99),], aes(x=UMAP_1, y=UMAP_2), color='gray', alpha=0.6) +
geom_point(data =  tmp[tmp$score >= quantile(tmp$score, probs=.99),], aes(x=UMAP_1, y=UMAP_2, color=nmf_clustering_11), alpha=1, size=3) + 
scale_color_manual(values=scanpy_colors_tmp) + theme_classic() + ggtitle('quantile99')

dev.off()

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/LSC_Scores_UMAP_IndGenes.pdf')
for (gene in c(positive_genes, negative_genes)){
	print(gene)
	tmp <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', gene))
	p <- ggplot(tmp, aes(x=UMAP_1, y= UMAP_2, color=tmp[[3]])) + 
				geom_point(size=1, alpha=0.6) +
				scale_color_distiller(palette = "Spectral", direction = -1, name = 'Score') + 
				theme_classic() + theme(legend.position="right", panel.background = element_blank(), 
				axis.ticks=element_blank(), axis.text=element_blank())
	if(gene %in% positive_genes){
		p <- p + ggtitle(paste0(gene, ' (positive)'))
	}else{
		p <- p + ggtitle(paste0(gene, ' (negative)'))
	}
	print(p)
}
dev.off()






get_umap_list <- function(gene_list, obj = combined_obj){
	data <- FetchData(obj, vars=c('UMAP_1', 'UMAP_2', gene_list))
	plot_list <- list()
	for (gene in gene_list){
		data <- data[order(data[[gene]]), ]
		p <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = .data[[gene]])) + 
		geom_point(size=1, alpha=0.9) + scale_color_distiller(palette="Greens", direction =1, name=gene) + 
		theme_classic() + theme(legend.position='bottom',  panel.background = element_blank(), 
				axis.ticks=element_blank(), axis.text=element_blank()) + 
		guides(color=guide_colorbar(barheight=0.2,barwidth = 2, label.position="bottom"))
		plot_list[[gene]] <- p
	}
	return(plot_list)
}

pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/LSC_signatures/LSC_signature_CONTROL.pdf')
DefaultAssay(combined_obj) <- "RNA"
tmp <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2', 'nCount_RNA'))
# tmp[tmp$nCount_RNA > 1000000, 'nCount_RNA'] <- 1000000
ggplot(tmp, aes(x=UMAP_1, y=UMAP_2, color=log(nCount_RNA))) + 
	geom_point(size=1, alpha=0.9) + scale_color_distiller(palette="Oranges", direction =1) + 
	theme_classic()

rdm_genes <-  sample(combined_obj@assays$RNA@var.features, 6)
a <- get_umap_list(rdm_genes)
cowplot::plot_grid(plotlist=a, nrow=2)
# FeaturePlot(combined_obj, features = c("ABI3" ,  "HOXB8" , "KRBOX1" ,"PTX3"  , "FOLR2" , "MCM7"), min.cutoff = "q9")
dev.off()



data <- FetchData(combined_obj, vars=c('UMAP_1', 'UMAP_2',"ABI3" ,  "HOXB8" , "KRBOX1" ,"PTX3"  , "FOLR2" , "MCM7"))

gene_list <- c("ABI3" ,  "HOXB8" , "KRBOX1" ,"PTX3"  , "FOLR2" , "MCM7")




library(Seurat)
library(ggplot2)
library(enrichR)


combined_obj <- readRDS( '/home/sevastopol/data/gserranos/SEACells_AML/Data/SEACells/AML_metacells_RPCA_integrated_ann.rds')

rank <- 11
res_clustering <- readRDS(paste0('/home/sevastopol/data/gserranos/SEACells_AML/Data/NMF/res_NMF_Rank_Integrated_Norm_',rank,'.rds'))
clustering_nmf <- setNames(as.data.frame(NMF::predict(res_clustering)), paste0('NMF_cluster_', rank))
combined_obj$nmf_clustering_11 <- clustering_nmf


C5_C2_pathways <- c(
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/c5.all.v2023.1.Hs.symbols.gmt'), 
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/c2.all.v2023.1.Hs.symbols.gmt'), 
	fgsea::gmtPathways('/home/sevastopol/data/gserranos/SEACells_AML/Data/Other/GOBP_RESPONSE_TO_DRUG.v2023.1.Hs.gmt')
	)



all_terms <- c(
	'GAL_LEUKEMIC_STEM_CELL_UP',
	'GAL_LEUKEMIC_STEM_CELL_DN',
	'GENTLES_LEUKEMIC_STEM_CELL_DN',
	'GENTLES_LEUKEMIC_STEM_CELL_UP',
	'GRAHAM_CML_QUIESCENT_VS_NORMAL_DIVIDING_UP',
	'GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP',
	'GRAHAM_CML_QUIESCENT_VS_CML_DIVIDING_UP')


DefaultAssay(combined_obj) <- 'RNA'

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



pdf('/home/sevastopol/data/gserranos/SEACells_AML/Plots/Signatures_quiescence.pdf')
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



