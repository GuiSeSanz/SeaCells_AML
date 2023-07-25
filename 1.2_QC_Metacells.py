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
# %matplotlib inline


sns.set_style("ticks")
matplotlib.rcParams["figure.figsize"] = [10, 10]
matplotlib.rcParams["figure.dpi"] = 100

# Load data
AML_SAMPLES = os.listdir("./Data/SEACells")


aml_dictionary = {
    "aml1": "SMD_AML_895521",
    "aml2": "AML34665",
    "aml3": "AML_CMM",
    "aml4": "AML_BZAA",
    "aml5": "AML37604",
    "aml6": "AML37742",
    "aml7": "AML38359",
    "aml8": "AML_43235",
    "aml9": "AML_43328",
    "aml10": "AML_29219",
    "aml11": "AML_34009",
}


aml_pattern = re.compile(r"(?<=_)[\w]+")

for aml_sample in AML_SAMPLES:
	print(aml_sample)
	DATA_PATH = "./Data/SEACells/" + aml_sample + "/"
	aml_data = sc.read(os.path.join(DATA_PATH, "AML_RNA_seacells.h5ad"))
	SEACell_ad = SEACells.core.summarize_by_SEACell(
		aml_data, SEACells_label="SEACell", summarize_layer="raw"
	)
	a = SEACell_ad.to_df()
	a.to_csv(f"{DATA_PATH}AML_seacells_RAW_data_colnames.csv", sep="\t", mode="w")




for aml_sample in AML_SAMPLES:
    print(aml_sample)
    PLOT_PATH = "./Plots/SEACells/" + aml_sample + "/"
    assert os.path.exists(PLOT_PATH), "Path does not exist"
    DATA_PATH = "./Data/SEACells/" + aml_sample + "/"
    assert os.path.exists(DATA_PATH), "Path does not exist"
    # Load data
    aml_data = sc.read(os.path.join(DATA_PATH, "AML_RNA_seacells.h5ad"))
    pd.DataFrame(aml_data.obsm["X_umap"]).set_index(aml_data.obs_names).to_csv(f'{DATA_PATH}/original_umap_coords.csv', sep="\t", mode="w")

    SEACells.plot.plot_2D(
        aml_data,
        key="X_umap",
        colour_metacells=False,
        save_as=f"{PLOT_PATH}X_umap_AML.pdf",
    )
    SEACells.plot.plot_2D(
        aml_data,
        key="X_umap",
        colour_metacells=True,
        save_as=f"{PLOT_PATH}X_umap_AML_colored.pdf",
    )
    SEACells.plot.plot_2D(
        aml_data,
        key="X_pca",
        colour_metacells=False,
        save_as=f"{PLOT_PATH}X_pca_AML.pdf",
    )
    SEACells.plot.plot_SEACell_sizes(
        aml_data, bins=10, save_as=f"{PLOT_PATH}SEACell_sizes_AML.pdf"
    )

    # SEACells.plot.plot_initialization(ad, model, save_as=f'{PLOT_PATH}Initialization_AML.pdf')

    aml_data.obs[["SEACell"]].to_csv(
        f"{DATA_PATH}AML_cells_2_SEACells.csv", sep="\t", mode="w"
    )
    df = pd.DataFrame.from_dict(
        Counter(aml_data.obs[["SEACell"]].loc[:, "SEACell"]),
        orient="index",
        columns=["SEACells"],
    )
    df.plot.bar()
    plt.savefig(f"{PLOT_PATH}Cell_distribution.pdf")
    umap = (
        pd.DataFrame(aml_data.obsm["X_umap"])
        .set_index(aml_data.obs_names)
        .join(aml_data.obs["SEACell"])
    )
    umap["SEACell"] = umap["SEACell"].astype("category")
    mcs = umap.groupby("SEACell").mean().reset_index()
    mcs.to_csv(f"{DATA_PATH}AML_seacells_coords_colnames.csv", sep="\t", mode="w")

    ### Metrics
    SEACell_purity = SEACells.evaluate.compute_celltype_purity(aml_data, "cluster")
    compactness = SEACells.evaluate.compactness(aml_data, "X_umap")
    separation = SEACells.evaluate.separation(aml_data, "X_umap", nth_nbr=1)

    plt.figure(figsize=(4, 4))
    sns.boxplot(data=SEACell_purity, y="cluster_purity")
    plt.title("cluster Purity")
    sns.despine()
    plt.savefig(f"{PLOT_PATH}AML_SEACell_purity.pdf", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(4, 4))
    sns.boxplot(data=compactness, y="compactness")
    plt.title("compactness")
    sns.despine()
    plt.savefig(f"{PLOT_PATH}AML_SEACell_compactness.pdf", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(4, 4))
    sns.boxplot(data=separation, y="separation")
    plt.title("cluster separation")
    sns.despine()
    plt.savefig(f"{PLOT_PATH}AML_SEACell_separation.pdf", bbox_inches="tight")
    plt.close()

    metrics = pd.concat(
        [compactness, separation, SEACell_purity], axis=1, ignore_index=True
    )
    metrics.columns = ["compactness", "separation", "cluster", "purity"]

    metrics.to_csv(f"{DATA_PATH}/AML_harmony_metacells_metrics.txt", sep=" ", mode="w")

    # # Save anndata
    aml_data.write(f"{DATA_PATH}AML_RNA_harmony_seacells.h5ad", compression="gzip")

    # # Summarize data
    SEACell_ad = SEACells.core.summarize_by_SEACell(
        aml_data, SEACells_label="SEACell", summarize_layer="raw"
    )
	a = SEACell_ad.to_df()
    a.to_csv(f"{DATA_PATH}AML_seacells_RAW_data_colnames.csv", sep="\t", mode="w")
    sc.pp.normalize_total(SEACell_ad)
    sc.pp.log1p(SEACell_ad)
    SEACell_ad.write(
        f"{DATA_PATH}/AML_RNA_harmony_seacells_metacells_normalized.h5ad",
        compression="gzip",
    )

    ### UMAP on SEACell matrix
    # SEACell_ad = sc.read(f'{DATA_PATH}/AML_RNA_harmony_seacells_metacells_normalized.h5ad')
    SEACell_ad.obs["# Single Cells"] = (
        aml_data.obs.groupby("SEACell").count().iloc[:, 0].loc[SEACell_ad.obs_names]
    )
    # SEACell_ad.obs['cluster'] = aml_data.obs.loc[SEACell_ad.obs_names, 'cluster']
    # SEACell_ad.obs['# Distinct Patients'] = aml_data.obs.groupby('SEACell').apply(lambda x: len(x['Patient'].unique())).loc[SEACell_ad.obs_names]
    SEACell_ad.obs["# Distinct Clusters"] = (
        aml_data.obs.groupby("SEACell")
        .apply(lambda x: len(x["cluster"].unique()))
        .loc[SEACell_ad.obs_names]
    )
    SEACell_ad.obs["compactness"] = metrics.loc[SEACell_ad.obs_names, "compactness"]
    SEACell_ad.obs["separation"] = metrics.loc[SEACell_ad.obs_names, "separation"]
    SEACell_ad.obs["purity"] = metrics.loc[SEACell_ad.obs_names, "purity"]

    SEACell_ad.obs.to_csv(f"{DATA_PATH}AML_seacells_metadata.csv", sep=" ")
    # SEACell_ad.X.to_csv(f'{DATA_PATH}AML_seacells_adata.csv', sep=' ')
    # np.savetxt(f'{DATA_PATH}AML_seacells_data.csv', SEACell_ad.X.todense(), delimiter=' ', header = SEACell_ad.var_names)

    a = SEACell_ad.to_df()
    a.to_csv(f"{DATA_PATH}AML_seacells_data_colnames.csv", sep="\t", mode="w")
    # SEACell_ad.obs.head()
    SEACell_ad.write(
        f"{DATA_PATH}AML_RNA_harmony_seacells_metacells_normalized_ann.h5ad",
        compression="gzip",
    )

    sc.pp.highly_variable_genes(SEACell_ad)
    sc.tl.pca(SEACell_ad, use_highly_variable=True)
    sc.pp.neighbors(SEACell_ad)
    sc.tl.umap(SEACell_ad)
    # sc.pl.umap(SEACell_ad, color='cluster', s=150)
    # sc.pl.umap(SEACell_ad, color='# Distinct Patients', s=150)
    sc.pl.umap(SEACell_ad, color="# Single Cells", s=150)
    plt.savefig(f"{PLOT_PATH}Single_cellsSC.pdf")
    plt.close()
    sc.pl.umap(SEACell_ad, color="compactness", s=150)
    plt.savefig(f"{PLOT_PATH}compactnessSC.pdf")
    plt.close()
    sc.pl.umap(SEACell_ad, color="separation", s=150)
    plt.savefig(f"{PLOT_PATH}separationSC.pdf")
    plt.close()
    sc.pl.umap(SEACell_ad, color="purity", s=150)
    plt.savefig(f"{PLOT_PATH}puritySC.pdf")
    plt.close()

    SEACell_ad.write(
        f"{DATA_PATH}aml_SCANPY_PROCESSED_seacells_metacells.h5ad", compression="gzip"
    )

    sc.pp.highly_variable_genes(SEACell_ad)
    plt.savefig(f"{PLOT_PATH}AML_SEACell_highly_variable_genes.pdf")
    plt.close()

    SEACell_ad.var["mt"] = SEACell_ad.var_names.str.startswith("MT-")
    SEACell_ad.var["rp"] = SEACell_ad.var_names.str.startswith("RPS")
    sc.pp.calculate_qc_metrics(
        SEACell_ad, qc_vars=["mt", "rp"], percent_top=None, log1p=False, inplace=True
    )
    sc.pl.violin(
        SEACell_ad,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rp"],
        jitter=0.4,
        multi_panel=True,
    )
    plt.savefig(f"{PLOT_PATH}AML_SEA_CellQC.pdf")
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(20, 6))
    sc.pl.scatter(
        SEACell_ad, x="total_counts", y="n_genes_by_counts", ax=ax[0], show=False
    )
    sc.pl.scatter(SEACell_ad, x="total_counts", y="pct_counts_mt", ax=ax[1], show=False)
    plt.savefig(f"{PLOT_PATH}AML_SEACell_n_genes_by_counts.pdf")
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(20, 6))
    sc.pp.highly_variable_genes(SEACell_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(SEACell_ad)
    plt.savefig(f"{PLOT_PATH}AML_SEACell_Highly_variable.pdf")
    plt.close()

    sc.tl.pca(SEACell_ad, svd_solver="arpack")
    fig, ax = plt.subplots(1, 2, figsize=(20, 6))
    sc.pl.pca(SEACell_ad, color="CST3")
    plt.savefig(f"{PLOT_PATH}AML_SEACell_PCA.pdf")
    plt.close()
    sc.pl.pca_variance_ratio(SEACell_ad, log=True)

    sc.pp.neighbors(SEACell_ad, n_neighbors=10, n_pcs=50)

    for res in [1, 1.2, 1.4, 1.6, 1.8, 2]:
        sc.tl.louvain(SEACell_ad, resolution=res, key_added="louvain_r" + str(res))
        sc.tl.leiden(SEACell_ad, resolution=res, key_added="leiden_r" + str(res))

    sc.tl.paga(SEACell_ad, groups="louvain_r1.6")
    sc.pl.paga(
        SEACell_ad, plot=True
    )  # remove `plot=False` if you want to see the coarse-grained graph
    plt.savefig(f"{PLOT_PATH}AML_SEACell_ClusterLouvain.pdf")
    plt.close()

    sc.tl.umap(SEACell_ad, init_pos="paga")
    sc.tl.umap(SEACell_ad)
    sc.pl.umap(SEACell_ad, color=["CD34", "HBB", "SATB1"])
    plt.savefig(f"{PLOT_PATH}AML_SEACell_UMAP_geneExpr.pdf")
    plt.close()

    sc.pl.umap(SEACell_ad, color=["louvain_r1.6"])
    plt.savefig(f"{PLOT_PATH}AML_SEACell_UMAP_leiden.pdf")
    plt.close()
    for res in [1, 1.2, 1.4, 1.6, 1.8, 2]:
        key = "louvain_r" + str(res)
        SEACell_ad.obs[key].to_csv(f"{DATA_PATH}AML_seacells_{key}.csv", sep="\t")
        key = "leiden_r" + str(res)
        SEACell_ad.obs[key].to_csv(f"{DATA_PATH}AML_seacells_{key}.csv", sep="\t")
