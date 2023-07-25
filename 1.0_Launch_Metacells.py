#! docker run -itd --name seacell_aml -v /home/sevastopol/data/gserranos/SEACells_AML:/AML/ alexthiery/seacells bash

############################
######    SEACells    ######
############################
# nohup docker exec seacell_aml bash /AML/1.0_Launch_Metacells.sh >/home/sevastopol/data/gserranos/SEACells_AML/AML_seacells_nohup.out 2>&1 &
### Perform SEACells directly on the 8 AML integrated patients data using Harmony reductions
import os
import re
import math
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from collections import Counter
# Some plotting aesthetics
#%matplotlib inline



def check_list_elements(list1, list2):
	common_elements = set(list1) & set(list2)
	unique_elements = set(list1) - set(list2)
	elements_in_both = len(common_elements)
	elements_only_in_list1 = len(unique_elements)
	elements_only_in_list2 = len(set(list2) - set(list1))
	print("Elements present in both lists:", elements_in_both)
	print("Elements only in list1:", elements_only_in_list1)
	print("Elements only in list2:", elements_only_in_list2)



sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [10, 10]
matplotlib.rcParams['figure.dpi'] = 100

# Load data  
files_aml = os.listdir('/AML/Data/')
files_aml = [file for file in files_aml if file.endswith('.h5ad')]


aml_dictionary = {
'aml1'  : 'SMD_AML_895521',
'aml2'  : 'AML34665',
'aml3'  : 'AML_CMM',
'aml4'  : 'AML_BZAA',
'aml5'  : 'AML37604',
'aml6'  : 'AML37742',
'aml7'  : 'AML38359',
'aml8'  : 'AML_43235',
'aml9'  : 'AML_43328',
'aml10' : 'AML_29219',
'aml11' : 'AML_34009'}



aml_pattern = re.compile(r'(?<=_)[\w]+')

for aml_sample in files_aml:
	sample_name = aml_pattern.search(aml_sample).group(0)
	print(aml_sample)
	PLOT_PATH = '/AML/Plots/SEACells/' + sample_name + '/'
	if not os.path.exists(PLOT_PATH):
		os.makedirs(PLOT_PATH)
	DATA_PATH = '/AML/Data/SEACells/'  +sample_name + '/'
	if not os.path.exists(DATA_PATH):
		os.makedirs(DATA_PATH)
	aml_sample_path = '/AML/Data/' + aml_sample
	assert os.path.exists(aml_sample_path), 'File not found.'
	ad = sc.read(aml_sample_path)
	# get the counts from the 10x h5 file
	raw_name = aml_dictionary.get(sample_name, 'Not Found')
	# ad2 = sc.read_10x_h5(f'/AML/Data/Raw_h5_objects/AML_matrices/{raw_name}_Count/outs/filtered_feature_bc_matrix.h5')
	# ad2.obs_names = [name.replace('-1', f'_{sample_name}') for name in ad2.obs_names]
	# set the new clinical classification
	# Secondary_samples = ['aml2', 'aml7', 'aml9', 'aml11']
	# category_2023 = ['Secondary' if sample in Secondary_samples  else 'Secondary' for sample in ad.obs['Patient']]
	# ad.obs['category_2023'] = pd.Series(category_2023, dtype="category", index=ad.obs_names)
	# Counter(ad.obs['Patient'])
	# Counter(ad.obs['category_2023'])
	# preprocessing
	# Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
	# This step should be performed after filtering 
	# check_list_elements(ad.obs_names.tolist(),ad2.obs_names.tolist())
	# check_list_elements(ad.var_names.tolist(),ad2.var_names.tolist())
	# ad2.var_names_make_unique()
	# raw_ad = sc.AnnData(ad2[ad.obs_names.tolist(), ad.var_names.tolist()].X)
	# raw_ad = sc.AnnData(ad.raw.X)
	# raw_ad.obs_names, raw_ad.var_names = ad.obs_names, ad.var_names
	# ad.raw = raw_ad
	# Normalize cells, log transform and compute highly variable genes
	sc.pp.normalize_per_cell(ad)
	sc.pp.log1p(ad)
	sc.pp.highly_variable_genes(ad, n_top_genes=3000)
	sc.tl.pca(ad, n_comps=50, use_highly_variable=True)
	sc.pp.neighbors(ad, n_pcs = 30, n_neighbors = 20)
	sc.tl.umap(ad)

	ad.obs['cluster'] = pd.Series(ad.obs['SCT_snn_res.0.5'], dtype="category")
	# Plot cell-types for reference
	plt.style.use('dark_background')
	sc.pl.scatter(ad, basis='umap', color='cluster', frameon=False, show=False)
	plt.savefig(f"{PLOT_PATH}UMAP_pyhton.pdf", bbox_inches="tight")


	### Running SEACells using harmony reduction
	# we recommended choosing one metacell for every 75 single-cells

	n_SEACells = math.ceil(ad.n_obs/75)
	build_kernel_on = 'X_umap' # key in ad.obsm to use for computing metacells
							# This would be replaced by 'X_svd' for ATAC data

	print(f'Creating model with {n_SEACells} metacells on {build_kernel_on}...')
	n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
	waypoint_proportion = 0.9 # Proportion of metacells to initialize using waypoint analysis, 
							# the remainder of cells are selected by greedy selection
	model = SEACells.core.SEACells(ad, 
					build_kernel_on=build_kernel_on, 
					n_SEACells=n_SEACells, 
					n_waypoint_eigs=n_waypoint_eigs,
					convergence_epsilon = 1e-5)

	print('Model created!')
	model.construct_kernel_matrix()
	# Initialize archetypes
	model.initialize_archetypes()

	# save model
	with open(f'{DATA_PATH}seacells_model_Untrained.obj', 'wb') as file_pi:
		pickle.dump(model, file_pi)

	print('Model training started!')
	model.fit(min_iter=10, max_iter=200)

	print('Model trained!')
	ad.__dict__['_raw'].__dict__['_var'] = ad.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
	# Save anndata
	ad.write(f'{DATA_PATH}AML_RNA_seacells.h5ad', compression='gzip')


	# # Save metacells metadata
	# txt_file = f'{DATA_PATH}AML_metacells.txt'
	# print('file saved to: ', txt_file)
	# ad.obs.to_csv(txt_file, sep=" ", mode="w")

	# # save model
	# with open(f'{DATA_PATH}AML_RNA_seacells_model_Trained.obj', 'wb') as file_pi:
	# 	pickle.dump(model, file_pi)

	# with open(f'{DATA_PATH}AML_RNA_seacells_model_Trained.obj', 'rb') as f:
	# 	model = pickle.load(f)

	# ad = sc.read(f'{DATA_PATH}AML_RNA_seacells.h5ad')


	# plt.style.use('dark_background')
	# model.plot_convergence(save_as=f'{PLOT_PATH}Convergence_AML.pdf')

	# SEACells.plot.plot_2D(ad, key='X_umap' ,        colour_metacells=False, save_as=f'{PLOT_PATH}X_umap_AMLy.pdf')
	# SEACells.plot.plot_2D(ad, key='X_umap' ,        colour_metacells=True,  save_as=f'{PLOT_PATH}X_umap_AML_colored.pdf')
	# # SEACells.plot.plot_2D(ad, key='X_harmony',      colour_metacells=False, save_as=f'{PLOT_PATH}X_harmony_AML.pdf')
	# # SEACells.plot.plot_2D(ad, key='X_non_int_UMAP', colour_metacells=False, save_as='/AML/Plots/SEACells/X_non_int_UMAP_AML.pdf')
	# SEACells.plot.plot_2D(ad, key='X_pca' ,         colour_metacells=False, save_as=f'{PLOT_PATH}X_pca_AML.pdf')
	# SEACells.plot.plot_SEACell_sizes(ad, bins=10, save_as=f'{PLOT_PATH}SEACell_sizes_AML.pdf')

	# SEACells.plot.plot_initialization(ad, model, save_as=f'{PLOT_PATH}Initialization_AML.pdf')

	# ad.obs[['SEACell']].to_csv(f'{DATA_PATH}AML_cells_2_SEACells.csv', sep='\t', mode='w')



	# umap = pd.DataFrame(ad.obsm['X_umap']).set_index(ad.obs_names).join(ad.obs['SEACell'])
	# umap['SEACell'] = umap['SEACell'].astype("category")
	# mcs = umap.groupby('SEACell').mean().reset_index()
	# mcs.to_csv(f'{DATA_PATH}AML_seacells_Integrated_coords_colnames.csv', sep='\t', mode='w')





# ### Metrics
# SEACell_purity  = SEACells.evaluate.compute_celltype_purity(ad, 'cluster')
# compactness = SEACells.evaluate.compactness(ad, 'X_harmony')
# separation = SEACells.evaluate.separation(ad, 'X_harmony',nth_nbr=1)

# plt.figure(figsize=(4,4))
# sns.boxplot(data=SEACell_purity, y='cluster_purity')
# plt.title('cluster Purity')
# sns.despine()
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_purity.pdf", bbox_inches="tight")
# plt.close()

# plt.figure(figsize=(4,4))
# sns.boxplot(data=compactness, y='compactness')
# plt.title('compactness')
# sns.despine()
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_compactness.pdf", bbox_inches="tight")
# plt.close()

# plt.figure(figsize=(4,4))
# sns.boxplot(data=separation, y='separation')
# plt.title('cluster separation')
# sns.despine()
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_separation.pdf", bbox_inches="tight")
# plt.close()

# metrics = pd.concat([compactness, separation, SEACell_purity], axis = 1, ignore_index = True)
# metrics.columns = ['compactness', 'separation', 'cluster', 'purity']

# metrics.to_csv('/AML/Data/SEACells/AML_harmony_metacells_metrics.txt', sep=" ", mode="w")




# plt.figure(figsize=(3,2))
# sns.distplot((model.A_.T > 0.1).sum(axis=1), kde=False)
# plt.title(f'Non-trivial (> 0.1) assignments per cell')
# plt.xlabel('# Non-trivial SEACell Assignments')
# plt.ylabel('# Cells')
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_Non-trivial_assignments.pdf", bbox_inches="tight")
# plt.close()


# plt.figure(figsize=(3,2))
# b = np.partition(model.A_.T, -5)    
# sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
# plt.title('Strength of top 5 strongest assignments')
# plt.xlabel('$n^{th}$ strongest assignment')
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_StringAssignment.pdf", bbox_inches="tight")
# plt.close()

# # # Save anndata
# ad.write('/AML/Data/SEACells/AML_RNA_harmony_seacells.h5ad', compression='gzip')

# # # Summarize data
# SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')
# sc.pp.normalize_total(SEACell_ad)
# sc.pp.log1p(SEACell_ad)


# SEACell_ad.write('/AML/Data/SEACells/AML_RNA_harmony_seacells_metacells_normalized2.h5ad', compression='gzip')


# ### UMAP on SEACell matrix

# SEACell_ad = sc.read('/AML/Data/SEACells/AML_RNA_harmony_seacells_metacells_normalized2.h5ad')
# SEACell_ad.obs['# Single Cells'] = ad.obs.groupby('SEACell').count().iloc[:,0].loc[SEACell_ad.obs_names]
# SEACell_ad.obs['cluster'] = ad.obs.loc[SEACell_ad.obs_names, 'cluster']
# SEACell_ad.obs['# Distinct Patients'] = ad.obs.groupby('SEACell').apply(lambda x: len(x['Patient'].unique())).loc[SEACell_ad.obs_names]
# SEACell_ad.obs['# Distinct Clusters'] = ad.obs.groupby('SEACell').apply(lambda x: len(x['cluster'].unique())).loc[SEACell_ad.obs_names]
# SEACell_ad.obs['compactness'] = metrics.loc[SEACell_ad.obs_names,'compactness']
# SEACell_ad.obs['separation'] = metrics.loc[SEACell_ad.obs_names,'separation']
# SEACell_ad.obs['purity'] = metrics.loc[SEACell_ad.obs_names,'purity']

# SEACell_ad.obs.to_csv('/AML/Data/SEACells/AML_seacells_metadata.csv', sep=' ')
# # SEACell_ad.X.to_csv('/AML/Data/SEACells/AML_seacells_metadata.csv', sep=' ')
# # np.savetxt('/AML/Data/SEACells/AML_seacells_data.csv', SEACell_ad.X.todense(), delimiter=' ', header = SEACell_ad.var_names)

# a = SEACell_ad.to_df()
# a.to_csv('/AML/Data/SEACells/AML_seacells_data_colnames.csv', sep='\t', mode='w')
# # SEACell_ad.obs.head()
# SEACell_ad.write('/AML/Data/SEACells/AML_RNA_harmony_seacells_metacells_normalized_ann.h5ad', compression='gzip')


# # sc.tl.pca(SEACell_ad, use_highly_variable=True)
# # sc.pp.neighbors(SEACell_ad, use_rep='X_pca')
# # sc.tl.umap(SEACell_ad)
# # sc.pl.umap(SEACell_ad, color='cluster', s=150)
# # sc.pl.umap(SEACell_ad, color='# Distinct Patients', s=150)
# # sc.pl.umap(SEACell_ad, color='# Single Cells', s=150)
# # sc.pl.umap(SEACell_ad, color='# Distinct Clusters', s=150)
# # sc.pl.umap(SEACell_ad, color='compactness', s=150)
# # sc.pl.umap(SEACell_ad, color='separation', s=150)
# # sc.pl.umap(SEACell_ad, color='purity', s=150)

# # SEACell_ad.write('./processed_data/aml_RNA_harmony_seacells_metacells.h5ad', compression='gzip')

# plt.figure(figsize=(3,5))
# sc.pl.highest_expr_genes(SEACell_ad, n_top=10, show=False)
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_highly_expressed_genes.pdf")
# plt.close()

# sc.pp.highly_variable_genes(SEACell_ad)
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_highly_variable_genes.pdf")
# plt.close()

# SEACell_ad.var['mt'] = SEACell_ad.var_names.str.startswith('MT-') 
# SEACell_ad.var['rp'] = SEACell_ad.var_names.str.startswith('RPS') 
# sc.pp.calculate_qc_metrics(SEACell_ad, qc_vars=['mt', 'rp'], percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(SEACell_ad, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rp'], jitter=0.4, multi_panel=True)
# plt.savefig("/AML/Plots/SEACells/AML_SEA_CellQC.pdf")
# plt.close()

# fig, ax = plt.subplots(1,2, figsize=(20,6))
# sc.pl.scatter(SEACell_ad, x='total_counts', y='n_genes_by_counts', ax=ax[0], show=False)
# sc.pl.scatter(SEACell_ad, x='total_counts', y='pct_counts_mt', ax=ax[1], show=False)
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_n_genes_by_counts.pdf")
# plt.close()



# fig, ax = plt.subplots(1,2, figsize=(20,6))
# sc.pp.highly_variable_genes(SEACell_ad, min_mean=0.0125, max_mean=3,  min_disp=0.5)
# sc.pl.highly_variable_genes(SEACell_ad)
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_Highly_variable.pdf")
# plt.close()

# sc.tl.pca(SEACell_ad, svd_solver='arpack')
# fig, ax = plt.subplots(1,2, figsize=(20,6))
# sc.pl.pca(SEACell_ad, color='CST3')
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_PCA.pdf")
# plt.close()
# sc.pl.pca_variance_ratio(SEACell_ad, log=True)

# sc.pp.neighbors(SEACell_ad, n_neighbors=10, n_pcs=50)

# for res in [1, 1.2, 1.4, 1.6, 1.8, 2]:
#     sc.tl.louvain(SEACell_ad, resolution=res, key_added='louvain_r'+str(res))
#     sc.tl.leiden(SEACell_ad, resolution=res, key_added='leiden_r'+str(res))

# sc.tl.paga(SEACell_ad, groups='louvain_r1.6')
# sc.pl.paga(SEACell_ad, plot=True)  # remove `plot=False` if you want to see the coarse-grained graph
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_graph.pdf")
# plt.close()

# sc.tl.umap(SEACell_ad, init_pos='paga')
# sc.tl.umap(SEACell_ad)
# sc.pl.umap(SEACell_ad, color=['CD34', 'HBB', 'SATB1'])
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_UMAP.pdf")
# plt.close()


# sc.pl.umap(SEACell_ad, color=['louvain_r1.6'])
# plt.savefig("/AML/Plots/SEACells/AML_SEACell_UMAP_leiden.pdf")
# plt.close()



# SEACell_ad.obs['leiden_r2'].to_csv('/AML/Data/SEACells/AML_seacells_LeidenClustering.csv', sep='\t')

# zeileis_26 = [
#     "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
#     "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
#     "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
#     "#f3e1eb", "#f6c4e1", "#f79cd4", "#7f7f7f", "#c7c7c7", "#1CE6FF", "#336600",
# ]

# godsnot_64 = [
#     # "#000000",  # remove the black, as often, we have black colored annotation
#     "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
#     "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
#     "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF",
#     "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
#     "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
#     "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
#     "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
#     "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
#     "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
#     "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
#     "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
#     "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
#     "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
#     "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
#     "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
# ]
