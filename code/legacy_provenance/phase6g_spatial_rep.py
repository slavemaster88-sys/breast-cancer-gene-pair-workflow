#!/usr/bin/env python
"""P1-4: replicate spatial ESV/sDDS gradients in 15 independent Visium BC
samples (GSE213688, fresh-frozen). Per sample: niche scores, ESV, sDDS;
boundary (malignant-low) vs core (malignant-high) comparison + Moran's I."""
import json
import joblib
import numpy as np
import pandas as pd
import scipy.stats as st
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
from anndata import read_h5ad
import squidpy as sq

D = "/workspace/SpaGI_BC/00_raw_data/Spatial_extra/extracted"
OUT = "/workspace/SpaGI_BC/06_clinical_validation"
P3 = "/workspace/SpaGI_BC/03_TME_extrinsic_SV"

sigs = json.load(open("/workspace/SpaGI_BC/05_spatial_niche/signature_genes.json"))
tri = pd.read_csv(f"{P3}/Tripartite_LR_SV_network.csv")
esv_module = sorted(set(tri["receptor"]) |
                    set(";".join(tri["receptor_upstream_SV_hubs"].dropna()).split(";")))
models = joblib.load("/workspace/SpaGI_BC/04_sDDS_model/sdds_models.joblib")

def score(adata, genes, name):
    g = [x for x in genes if x in adata.var_names]
    if not g:
        adata.obs[name] = 0.0
        return
    sc.tl.score_genes(adata, g, score_name=name, use_raw=False)

def sdds_spots(adata, drug):
    mdl = models[drug]
    Z = adata.raw.X if adata.raw is not None else adata.X
    import scipy.sparse as sp
    Xd = Z.toarray() if sp.issparse(Z) else np.asarray(Z)
    genes = adata.raw.var_names if adata.raw is not None else adata.var_names
    Zdf = pd.DataFrame(Xd, index=adata.obs_names, columns=genes)
    mu, sd = Zdf.mean(), Zdf.std() + 1e-9
    Zz = (Zdf - mu) / sd
    feats = pd.DataFrame(index=adata.obs_names)
    for c in mdl["features"]:
        a, b = c.split("_")[0], c.split("_")[1]
        if a in Zz.columns and b in Zz.columns:
            feats[c] = (Zz[a] + Zz[b]) / 2
    F = pd.DataFrame(0.0, index=adata.obs_names, columns=mdl["features"])
    cf = [f for f in mdl["features"] if f in feats.columns]
    if cf:
        F[cf] = feats[cf]
    return pd.Series(mdl["xgboost"].predict_proba(F)[:, 1], index=adata.obs_names)

import glob, os
h5s = sorted(glob.glob(f"{D}/*_filtered_feature_bc_matrix.h5"))
rows = []
for h5 in h5s:
    sample = os.path.basename(h5).split("_")[1]
    try:
        adata = sc.read_visium(os.path.dirname(h5), prefix=f"{os.path.basename(h5).split('_')[0]}_{sample}_")
    except Exception:
        # fallback: build from h5 + positions
        adata = sc.read_10x_h5(h5)
        adata.var_names_make_unique()
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    for nm, key in [("mal", "Malignant"), ("myCAF", "myCAF"), ("SPP1_TAM", "SPP1_TAM"),
                    ("CD8_T", "CD8_T")]:
        score(adata, sigs.get(key, []), f"s_{nm}")
    score(adata, esv_module, "esv_mod")
    tme_w = (adata.obs["s_myCAF"] + adata.obs["s_SPP1_TAM"]) / 2
    tw = (tme_w - tme_w.min()) / (tme_w.max() - tme_w.min() + 1e-9)
    adata.obs["ESV"] = adata.obs["esv_mod"] * tw
    adata.obs["sDDS_Paclitaxel"] = sdds_spots(adata, "Paclitaxel").values
    # core vs boundary by malignant score quartiles
    q1, q3 = adata.obs["s_mal"].quantile([0.25, 0.75])
    core = adata.obs["ESV"][adata.obs["s_mal"] >= q3]
    bnd = adata.obs["ESV"][adata.obs["s_mal"] <= q1]
    core_p = adata.obs["sDDS_Paclitaxel"][adata.obs["s_mal"] >= q3]
    bnd_p = adata.obs["sDDS_Paclitaxel"][adata.obs["s_mal"] <= q1]
    # Moran's I for ESV
    try:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
        sq.gr.spatial_autocorr(adata, obsm_keys=False, genes=["ESV"], mode="moran")
        mi = float(adata.uns["moranI"]["ESV"]["I"]) if "ESV" in adata.uns["moranI"].index else np.nan
    except Exception as e:
        mi = np.nan
    rows.append({"sample": sample,
                 "n_spots": adata.n_obs,
                 "ESV_core": core.mean(), "ESV_boundary": bnd.mean(),
                 "sDDS_core": core_p.mean(), "sDDS_boundary": bnd_p.mean(),
                 "moransI_ESV": mi})
    print(sample, f"n={adata.n_obs} ESV core/bnd={core.mean():.2f}/{bnd.mean():.2f} "
          f"sDDS core/bnd={core_p.mean():.3f}/{bnd_p.mean():.3f} MoranI={mi:.2f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/spatial_replication_gse213688.csv", index=False)
df_l = df.dropna(subset=["ESV_core"])
u, p_esv = st.wilcoxon(df_l.ESV_boundary, df_l.ESV_core)
u2, p_sdds = st.wilcoxon(df_l.sDDS_boundary, df_l.sDDS_core)
print(f"\nPaired Wilcoxon (n={len(df_l)} samples): "
      f"ESV boundary vs core P={p_esv:.2e}; sDDS_Paclitaxel boundary vs core P={p_sdds:.2e}")
print("P1_4_DONE")
