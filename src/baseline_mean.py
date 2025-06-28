"""
memory-safe mean-expression baseline
------------------------------------
Produces pred_validation_meanexp.h5ad without loading the full training
matrix into RAM (works on 8 GB laptops).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ─── paths ────────────────────────────────────────────────────────────────
DATA_DIR   = Path("data")
TRAIN      = DATA_DIR / "adata_Training.h5ad"
VAL_ROSTER = DATA_DIR / "pert_counts_Validation.csv"
GENES_CSV  = DATA_DIR / "gene_names.csv"

# ─── params ───────────────────────────────────────────────────────────────
PERT_COL = "target_gene"          # column in .obs and in the roster CSV
CTRL_KEYWORDS = ("ctrl", "ntc")   # how control rows are labelled (case-insensitive)

# ─── load roster & gene list ──────────────────────────────────────────────
val_perts = pd.read_csv(VAL_ROSTER)[PERT_COL].tolist()
gene_list = pd.read_csv(GENES_CSV, header=None)[0].tolist()
n_genes   = len(gene_list)

# ─── open training AnnData lazily ─────────────────────────────────────────
adata = sc.read_h5ad(TRAIN, backed="r")
obs_col = adata.obs[PERT_COL]

unique_perts = obs_col.unique()
print(f"{len(unique_perts)} unique perturbations found in training set")

# pre-allocate mean-matrix (perturbations x genes)
mean_mat = {}

# ─── iterate perturbations one by one ─────────────────────────────────────
for pert in unique_perts:
    mask = obs_col == pert
    sub  = adata[mask]                   # backed slice (~1 000 cells)
    # sub.X is CSR → mean(axis=0) is cheap
    mean_vec = np.asarray(sub.X.mean(axis=0)).ravel()   # 1-D float64
    mean_mat[pert] = mean_vec.astype("float32")
    # slice is freed when 'sub' goes out of scope

print("Finished computing per-perturbation means")

# ─── fallback vector (first control match or global mean) ─────────────────
ctrl_candidates = [p for p in mean_mat if any(k in p.lower() for k in CTRL_KEYWORDS)]
if ctrl_candidates:
    ctrl_vec = mean_mat[ctrl_candidates[0]]
else:
    # global mean (rarely needed)
    ctrl_vec = np.mean(list(mean_mat.values()), axis=0).astype("float32")

# ─── build prediction matrix for validation roster ───────────────────────
rows = [mean_mat.get(p, ctrl_vec) for p in val_perts]
pred_matrix = np.vstack(rows)

pred = ad.AnnData(
    X   = pred_matrix,
    obs = pd.DataFrame({PERT_COL: val_perts}),
    var = pd.DataFrame(index=gene_list),
    uns = {"method": "mean_expression_stream"}
)

OUT = Path("pred_validation_meanexp.h5ad")
pred.write_h5ad(OUT, compression="gzip")
print("Saved", OUT, "shape=", pred.shape)
