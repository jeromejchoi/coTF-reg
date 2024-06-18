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

#%%
mc.ut.logger().setLevel(logging.DEBUG)
# %%
# shutil.rmtree("/ua/choi267/Oligo project/GE_data/SEA-AD/MTG/full/full/metacells", ignore_errors=True)
# shutil.rmtree("/ua/choi267/Oligo project/GE_data/SEA-AD/MTG/full/full", ignore_errors=True)
os.makedirs("/figures", exist_ok=True)
os.makedirs("/ROSMAP/final", exist_ok=True)

# %%
full = ad.read_h5ad("rosmap_oligo_control.h5ad")
full.X = full.X.astype(np.float32)
full.X.sum_duplicates()
mc.ut.top_level(full)

mc.ut.set_name(full, "rosmap_oligo_control")
print(f"Full: {full.n_obs} cells, {full.n_vars} genes")


protein_gene = pd.read_csv("protein_coding_genes.csv")
# full_df_feature = pd.DataFrame(anndata_CTL.var)
full_df_feature = pd.DataFrame(full.var.index)
full_df_feature.columns = ['genes']
full_df_feature_df = full_df_feature[full_df_feature['genes'].isin(protein_gene['gene'])]
full_df_feature_s = full_df_feature.squeeze()
protein_gene_mask = full_df_feature_s.isin(protein_gene['gene'])

# # %%
full_df = full[:,full.var.index.isin(full_df_feature_df['genes'])]
full_df.X.sum_duplicates()

# %%
PROPERLY_SAMPLED_MIN_CELL_TOTAL = 800
PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000

# %%
total_umis_per_cell = mc.ut.get_o_numpy(full_df, "__x__", sum=True)
plot = sb.displot(total_umis_per_cell, log_scale=(10, None))
plot.set(xlabel="UMIs", ylabel="Density", yticks=[])

plot.refline(x=PROPERLY_SAMPLED_MIN_CELL_TOTAL, color="darkgreen")
plot.refline(x=PROPERLY_SAMPLED_MAX_CELL_TOTAL, color="crimson")

plt.savefig("rosmap_oligo_control.svg")

too_small_cells_count = np.sum(total_umis_per_cell < PROPERLY_SAMPLED_MIN_CELL_TOTAL)
too_large_cells_count = np.sum(total_umis_per_cell > PROPERLY_SAMPLED_MAX_CELL_TOTAL)

total_umis_per_cell = mc.ut.get_o_numpy(full_df, name="__x__", sum=True)
too_small_cells_percent = 100.0 * too_small_cells_count / full_df.n_obs
too_large_cells_percent = 100.0 * too_large_cells_count / full_df.n_obs

print(
    f"Will exclude {too_small_cells_count} ({too_small_cells_percent:.2f}%%) cells"
    f" with less than {PROPERLY_SAMPLED_MIN_CELL_TOTAL} UMIs"
)
print(
    f"Will exclude {too_large_cells_count} ({too_large_cells_percent:.2f}%%) cells"
    f" with more than {PROPERLY_SAMPLED_MAX_CELL_TOTAL} UMIs"
)


# %%
plt.boxplot(total_umis_per_cell)

# %%
min(total_umis_per_cell)
max(total_umis_per_cell)
np.mean(total_umis_per_cell)
np.median(total_umis_per_cell)

# %%
EXCLUDED_GENE_NAMES = []
EXCLUDED_GENE_PATTERNS = ["MT-.*"]  # Mytochondrial.

# %%
import scipy as sp


# %%
mc.pl.exclude_genes(
    full_df,
    excluded_gene_names=EXCLUDED_GENE_NAMES, 
    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS,
    random_seed=123456,
)

# %%
mc.tl.compute_excluded_gene_umis(full_df)

# %%
PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25

# %%
excluded_umis_fraction_regularization = 1e-3  # Avoid 0 values in log scale plot.
excluded_umis_per_cell = mc.ut.get_o_numpy(full_df, "excluded_umis")
excluded_umis_fraction_per_cell = excluded_umis_per_cell / total_umis_per_cell

excluded_umis_fraction_per_cell += excluded_umis_fraction_regularization
plot = sb.displot(excluded_umis_fraction_per_cell, log_scale=(10, None))
excluded_umis_fraction_per_cell -= excluded_umis_fraction_regularization

plot.set(xlabel="Fraction of excluded gene UMIs", ylabel="Density", yticks=[])
plot.refline(x=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION, color="crimson")

plt.savefig("rosmap_oligo_control_cell_excluded_umis_fraction.svg")

too_excluded_cells_count = np.sum(
    excluded_umis_fraction_per_cell > PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
)
too_excluded_cells_fraction = too_excluded_cells_count / full_df.n_obs

print(
    f"Will exclude {too_excluded_cells_count} ({100 * too_excluded_cells_fraction:.2f}%) cells"
    f" with more than {100 * PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION:.2f}% excluded gene UMIs"
)

# %%
mc.pl.exclude_cells(
    full_df,
    properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
    properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
    properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
)

# %%
clean = mc.pl.extract_clean_data(full_df, name="rosmap_oligo_control.clean")
mc.ut.top_level(clean)
print(f"Clean: {clean.n_obs} cells, {clean.n_vars} genes")

# %%
cells = clean
clean = None  # Allow it to be gc-ed
mc.ut.set_name(cells, "rosmap_oligo_control.cells")
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
mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456)

# %%
metacells = \
    mc.pl.collect_metacells(cells, name="rosmap_oligo_control.preliminary.metacells", random_seed=123456)
print(f"Preliminary: {metacells.n_obs} metacells, {metacells.n_vars} genes")

# %%
# mc.pl.compute_for_mcview(adata=cells, gdata=metacells, random_seed=123456)

# %%
# Assign a single value for each metacell based on the cells.
mc.tl.convey_obs_to_group(
    adata=cells, gdata=metacells,
    property_name="subject", to_property_name="subject",
    method=mc.ut.most_frequent  # This is the default, for categorical data
)

# Compute the fraction of cells with each possible value in each metacell:
mc.tl.convey_obs_fractions_to_group(
    adata=cells, gdata=metacells,
    property_name="subject", to_property_name="subject"
)


# %%
# type_of_metacell = np.array(metacell_types_csv["cell_type"])
# mc.ut.set_o_data(metacells, "type", type_of_metacell)

# extended_type_of_metacell = pd.Series(
#     list(type_of_metacell) + ["Outliers"],
#     index=list(metacell_types_csv["metacell"]) + ["Outliers"]
# )

# metacell_of_cell = cells.obs["metacell_name"]
# type_of_cell = np.array(extended_type_of_metacell[metacell_of_cell])
# mc.ut.set_o_data(cells, "type", type_of_cell)

# %%
cells.write_h5ad("rosmap_oligo_control_cells.h5ad")
metacells.write_h5ad("rosmap_oligo_control_metacells.h5ad")


