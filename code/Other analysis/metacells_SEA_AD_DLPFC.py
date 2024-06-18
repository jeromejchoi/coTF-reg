# %%
import anndata as ad             # For reading/writing AnnData files
import matplotlib.pyplot as plt  # For plotting
import numba
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sb             # For plotting
import scipy.sparse as sp        # For sparse matrices
import shutil                    # for filesystem operations
from math import hypot           # For plotting
import logging

mc.ut.logger().setLevel(logging.DEBUG)
# %%
full = ad.read_h5ad("SEA_AD_DLPFC.h5ad")
full.X.sum_duplicates()
# full_CTL = full[full.obs["disease"] == "normal", :]
status_mask = full.obs["disease"] == "normal"
full_CTL = mc.ut.slice(full, obs=status_mask)
mc.ut.top_level(full_CTL)
mc.ut.set_name(full_CTL, "SEA_AD_oligo_DLPFC")
print(f"Full: {full_CTL.n_obs} cells, {full_CTL.n_vars} genes")

# %%
print(full_CTL.X.__class__)

# %%
protein_gene = pd.read_csv("protein_coding_genes.csv") # Your path here

# %%
full_CTL_feature = pd.DataFrame(full_CTL.var["feature_name"])


# %%
full_CTL_feature_df = full_CTL_feature[full_CTL_feature.feature_name.isin(protein_gene['gene'])]


# %%
full_CTL_feature_s = full_CTL_feature.squeeze()

# %%
protein_gene_mask = full_CTL_feature_s.isin(protein_gene['gene'])

# %%
full_CTL_df = full_CTL[:,full_CTL.var['feature_name'].isin(list(full_CTL_feature_df['feature_name']))]
# full_CTL_df = mc.ut.slice(full_CTL, vars=protein_gene_mask)

# %%
# full_CTL_df.obs["value"] = 0


# %%
PROPERLY_SAMPLED_MIN_CELL_TOTAL = 800
PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000

# %%
total_umis_per_cell = mc.ut.get_o_numpy(full_CTL_df, "__x__", sum=True)
plot = sb.displot(total_umis_per_cell, log_scale=(10, None))
plot.set(xlabel="UMIs", ylabel="Density", yticks=[])

plot.refline(x=PROPERLY_SAMPLED_MIN_CELL_TOTAL, color="darkgreen")
plot.refline(x=PROPERLY_SAMPLED_MAX_CELL_TOTAL, color="crimson")

plt.savefig("cell_total_umis.svg")

too_small_cells_count = np.sum(total_umis_per_cell < PROPERLY_SAMPLED_MIN_CELL_TOTAL)
too_large_cells_count = np.sum(total_umis_per_cell > PROPERLY_SAMPLED_MAX_CELL_TOTAL)

total_umis_per_cell = mc.ut.get_o_numpy(full_CTL_df, name="__x__", sum=True)
too_small_cells_percent = 100.0 * too_small_cells_count / full_CTL_df.n_obs
too_large_cells_percent = 100.0 * too_large_cells_count / full_CTL_df.n_obs

print(
    f"Will exclude {too_small_cells_count} ({too_small_cells_percent:.2f}%%) cells"
    f" with less than {PROPERLY_SAMPLED_MIN_CELL_TOTAL} UMIs"
)
print(
    f"Will exclude {too_large_cells_count} ({too_large_cells_percent:.2f}%%) cells"
    f" with more than {PROPERLY_SAMPLED_MAX_CELL_TOTAL} UMIs"
)

# %%
# pip show metacells

# %%
EXCLUDED_GENE_NAMES = []
EXCLUDED_GENE_PATTERNS = ["MT-.*"]  # Mytochondrial.

# %%
mc.pl.exclude_genes(
    full_CTL_df,
    excluded_gene_names=EXCLUDED_GENE_NAMES, 
    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS,
    random_seed=123456,
)

# %%
mc.tl.compute_excluded_gene_umis(full_CTL_df)

# %%
PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25

# %%
excluded_umis_fraction_regularization = 1e-3  # Avoid 0 values in log scale plot.
excluded_umis_per_cell = mc.ut.get_o_numpy(full_CTL_df, "excluded_umis")
excluded_umis_fraction_per_cell = excluded_umis_per_cell / total_umis_per_cell

excluded_umis_fraction_per_cell += excluded_umis_fraction_regularization
plot = sb.displot(excluded_umis_fraction_per_cell, log_scale=(10, None))
excluded_umis_fraction_per_cell -= excluded_umis_fraction_regularization

plot.set(xlabel="Fraction of excluded gene UMIs", ylabel="Density", yticks=[])
plot.refline(x=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION, color="crimson")

plt.savefig("cell_excluded_umis_fraction.svg")

too_excluded_cells_count = np.sum(
    excluded_umis_fraction_per_cell > PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
)
too_excluded_cells_fraction = too_excluded_cells_count / full_CTL_df.n_obs

print(
    f"Will exclude {too_excluded_cells_count} ({100 * too_excluded_cells_fraction:.2f}%) cells"
    f" with more than {100 * PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION:.2f}% excluded gene UMIs"
)

# %%
mc.pl.exclude_cells(
    full_CTL_df,
    properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
    properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
    properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
)

# %%
clean = mc.pl.extract_clean_data(full_CTL_df, name="SEA_AD_DLPFC.clean")
mc.ut.top_level(clean)
print(f"Clean: {clean.n_obs} cells, {clean.n_vars} genes")

# %%
cells = clean
clean = None  # Allow it to be gc-ed
mc.ut.set_name(cells, "SEA_AD_DLPFC.preliminary.cells")
print(f"Input: {cells.n_obs} cells, {cells.n_vars} genes")

# %%
LATERAL_GENE_NAMES = []
LATERAL_GENE_PATTERNS = ["RP[LS].*"]  # Ribosomal

# %%
# This will mark as "lateral_gene" any genes that match the above, if they exist in the clean dataset.
mc.pl.mark_lateral_genes(
    cells,
    lateral_gene_names=LATERAL_GENE_NAMES,
    lateral_gene_patterns=LATERAL_GENE_PATTERNS,
)

lateral_gene_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
lateral_gene_names = set(cells.var_names[lateral_gene_mask])
print(sorted([
    name for name in lateral_gene_names
    if not name.startswith("RPL") and not name.startswith("RPS")
]))
print(f"""and {len([
    name for name in lateral_gene_names if name.startswith("RPL") or name.startswith("RPS")
])} RP[LS].* genes""")

# %%
NOISY_GENE_NAMES = []
mc.pl.mark_noisy_genes(cells, noisy_gene_names=NOISY_GENE_NAMES)

# %%
max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
# Or, if running out of memory manually override:
# max_paralle_piles = ...
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)

# %%
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456)

# %%
metacells = \
    mc.pl.collect_metacells(cells, name="SEA_AD_DLPFC.preliminary.metacells", random_seed=123456)
print(f"Preliminary: {metacells.n_obs} metacells, {metacells.n_vars} genes")

# %%
with mc.ut.progress_bar():
    mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=123456)

# %%
min_long_edge_size = 4
umap_x = mc.ut.get_o_numpy(metacells, "x")
umap_y = mc.ut.get_o_numpy(metacells, "y")
umap_edges = sp.coo_matrix(mc.ut.get_oo_proper(metacells, "obs_outgoing_weights"))
sb.set()
plot = sb.scatterplot(x=umap_x, y=umap_y, s=10)
for (
    source_index, target_index, weight
) in zip(
    umap_edges.row, umap_edges.col, umap_edges.data
):
    source_x = umap_x[source_index]
    target_x = umap_x[target_index]
    source_y = umap_y[source_index]
    target_y = umap_y[target_index]
    if hypot(target_x - source_x, target_y - source_y) >= min_long_edge_size:
        plt.plot([source_x, target_x], [source_y, target_y],
                 linewidth=weight * 2, color='indigo')
plt.show()

# %%
metacells.write_h5ad("SEA_AD_DLPFC_metacells.h5ad")
cells.write_h5ad("SEA_AD_DLPFC_cells.h5ad")

