# add_non_targeting.py
import anndata as ad, pandas as pd, numpy as np
IN  = "pred_validation_meanexp.h5ad"
OUT = "pred_validation_meanexp_ntc.h5ad"
PERT_COL = "target_gene"

adata = ad.read_h5ad(IN)
ctrl_vec = np.asarray(adata.X.mean(axis=0)).astype("float32")

ctrl = ad.AnnData(
    X   = ctrl_vec[np.newaxis, :],
    obs = pd.DataFrame({PERT_COL: ["non-targeting"]},
                       index=["non-targeting"]),
    var = adata.var.copy(),
)
adata2 = ad.concat([adata, ctrl], axis=0)
adata2.write_h5ad(OUT, compression="gzip")
print("Wrote", OUT, "shape", adata2.shape)
