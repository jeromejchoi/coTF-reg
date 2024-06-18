import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import gc
from scipy import sparse
import logging
import time

# Function to save AnnData object
def save_checkpoint(adata, filename):
    adata.write(filename)
    logging.info(f"Checkpoint saved: {filename}")
    
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

# Set figure parameters
sc.settings.set_figure_params(dpi=80)

# Load the data
logging.info("Loading data...")
adata = ad.read_h5ad("SEA_AD_DLPFC.h5ad")

# Filter for normal cells
logging.info("Filtering for normal cells...")
status_mask = adata.obs["disease"] == "normal"
full_CTL = adata[status_mask]

del adata
gc.collect()

# Load protein coding genes
logging.info("Loading protein coding genes...")
protein_gene = pd.read_csv("protein_coding_genes.csv")
full_CTL_feature = pd.DataFrame(full_CTL.var["feature_name"])
full_CTL_feature_df = full_CTL_feature[full_CTL_feature.feature_name.isin(protein_gene['gene'])]
full_CTL_df = full_CTL[:, full_CTL.var['feature_name'].isin(list(full_CTL_feature_df['feature_name']))]

# # Convert matrix to sparse format
# logging.info("Converting matrix to sparse format...")
# full_CTL_df.X = sparse.csr_matrix(full_CTL_df.X)

# Normalize and logarithmize
logging.info("Normalizing and logarithmizing...")
sc.pp.normalize_per_cell(full_CTL_df, counts_per_cell_after=1e4)
sc.pp.log1p(full_CTL_df)

# Store raw counts
full_CTL_df.raw = full_CTL_df

# Compute variable genes
logging.info("Computing variable genes...")
full_CTL_df.var['gene_ids'] = full_CTL_df.var['feature_name']
sc.pp.highly_variable_genes(full_CTL_df, min_mean=0.0125, max_mean=3, min_disp=0.5)
logging.info(f"Highly variable genes: {sum(full_CTL_df.var.highly_variable)}")

# Plot variable genes
sc.pl.highly_variable_genes(full_CTL_df)

# Subset for variable genes
logging.info("Subsetting for variable genes...")
full_CTL_df = full_CTL_df[:, full_CTL_df.var['highly_variable']]

# Scale data
logging.info("Scaling data...")
sc.pp.scale(full_CTL_df, max_value=10)

# PCA
logging.info("Running PCA...")
sc.tl.pca(full_CTL_df, svd_solver='arpack')
sc.pl.pca(full_CTL_df, color='Subclass', components = ['1,2','3,4','5,6','7,8'], ncols=2)
sc.pl.pca_loadings(full_CTL_df, components=[1,2,3,4,5,6,7,8])
sc.pl.pca_variance_ratio(full_CTL_df, log=True, n_pcs = 50)

# Compute neighbors
logging.info("Computing neighbors...")
sc.pp.neighbors(full_CTL_df, n_pcs = 30, n_neighbors = 20)

# UMAP
logging.info("Computing UMAP...")
sc.tl.umap(full_CTL_df)

logging.info("Script completed successfully.")

# Save the AnnData object with UMAP coordinates
save_checkpoint(full_CTL_df, 'full_CTL_df_with_umap.h5ad')
