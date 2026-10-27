#!/usr/bin/env python
"""P0-2 step 2: NNLS deconvolution of TCGA/METABRIC bulk RNA-seq using the
Phase-1 scRNA reference; recompute ESV with deconvolved myCAF/SPP1-TAM
fractions; rerun Cox models; compare with marker-proxy ESV."""
import json
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
from scipy.optimize import nnls
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

RAW = "/workspace/SpaGI_BC/00_raw_data/Clinical"
P1 = "/workspace/SpaGI_BC/01_scRNA_processed"
P3 = "/workspace/SpaGI_BC/03_TME_extrinsic_SV"
OUT = "/workspace/SpaGI_BC/06_clinical_validation"

tri = pd.read_csv(f"{P3}/Tripartite_LR_SV_network.csv")
esv_module = sorted(set(tri["receptor"]) |
                    set(";".join(tri["receptor_upstream_SV_hubs"].dropna()).split(";")))

# ---------- scRNA reference profiles ----------
a = sc.read_h5ad(f"{P1}/BC_scRNA_annotated.h5ad", backed="r")
hv = [g for g in a.var[a.var["highly_variable"]].index
      if not g.startswith(("MT-", "RPS", "RPL"))][:1000]
basis_genes = sorted(set(hv) | set(esv_module) |
                     {"COL1A1", "DCN", "FAP", "TAGLN", "SPP1", "CD68", "LYZ",
                      "CD3D", "CD8A", "GZMB", "ESR1", "ERBB2"})
basis_genes = [g for g in basis_genes if g in a.var_names]
print("basis genes in scRNA:", len(basis_genes), flush=True)

pos = a.var_names.get_indexer(basis_genes)
order = np.argsort(pos)          # h5py fancy indexing needs increasing order
inv = np.empty_like(order)
inv[order] = np.arange(len(order))
X = np.empty((a.n_obs, len(basis_genes)), dtype=np.float32)
CH = 20000
for s in range(0, a.n_obs, CH):
    e = min(s + CH, a.n_obs)
    blk = np.asarray(a.X[s:e, pos[order]])
    X[s:e][:, inv] = np.log1p(blk)
obs = a.obs
ct = obs["w_celltype_minor"].astype(str).values
spp1 = pd.to_numeric(obs["tme_SPP1_TAM"], errors="coerce").values
# refine SPP1 TAM: macrophage/monocyte cells with top-quartile SPP1 score
mac = np.isin(ct, ["Macrophage", "Monocyte", "Cycling_Myeloid"])
thr = np.nanquantile(spp1[mac], 0.75)
lab = ct.copy()
lab[mac & (spp1 >= thr)] = "SPP1_TAM"
lab[mac & (spp1 < thr)] = "Macrophage_other"
# collapse to deconvolution basis types
collapse = {
    "Cancer LumA SC": "Malignant", "Cancer LumB SC": "Malignant",
    "Cancer Basal SC": "Malignant", "Cancer Her2 SC": "Malignant",
    "Cancer Cycling": "Malignant", "Myoepithelial": "Malignant",
    "Luminal Progenitors": "Normal_Epi", "Mature Luminal": "Normal_Epi",
    "T cells CD4+": "CD4_T", "T cells CD8+": "CD8_T", "Cycling T-cells": "CD8_T",
    "NKT cells": "CD8_T", "NK cells": "NK",
    "B cells Memory": "B_cells", "B cells Naive": "B_cells",
    "Plasmablasts": "Plasma",
    "CAFs myCAF-like": "myCAF", "CAFs MSC iCAF-like": "iCAF",
    "Endothelial ACKR1": "Endo", "Endothelial CXCL12": "Endo",
    "Endothelial RGS5": "Endo", "Endothelial Lymphatic LYVE1": "Endo",
    "PVL Differentiated": "PVL", "PVL Immature": "PVL",
    "DCs": "DC",
}
lab = pd.Series(lab).replace(collapse).values
types = sorted(set(lab))
print("reference types:", types, flush=True)
ref = np.stack([X[lab == t].mean(axis=0) for t in types], axis=1)  # genes x types
ref_df = pd.DataFrame(ref, index=basis_genes, columns=types)
ref_df.to_csv(f"{OUT}/deconv_reference_profiles.csv")

# scale: z-score genes across cell types (CIBERSORT-style)
ref_z = (ref - ref.mean(axis=1, keepdims=True)) / (ref.std(axis=1, keepdims=True) + 1e-9)

def deconvolute(bulk):
    """bulk: genes x samples (log2). Returns fractions DataFrame."""
    bulk = bulk.dropna(axis=0, how="any")
    genes = [g for g in basis_genes if g in bulk.index]
    B = bulk.loc[genes].T.values.astype(float)  # samples x genes
    R = ref_z[[basis_genes.index(g) for g in genes]]  # genes x types
    # z-score bulk genes across samples
    Bz = (B - B.mean(0)) / (B.std(0) + 1e-9)
    ok = np.isfinite(R).all(axis=1) & np.isfinite(Bz).all(axis=0)
    R, Bz = R[ok], Bz[:, ok]
    fracs = np.zeros((Bz.shape[0], R.shape[1]))
    for i in range(Bz.shape[0]):
        coef, _ = nnls(R, Bz[i])
        fracs[i] = coef / (coef.sum() + 1e-9)
    return pd.DataFrame(fracs, index=bulk.columns, columns=types)

results = {}
for tag, age_col, stage_kind, sub_col in [
        ("TCGA", "AGE", "tcga", None),
        ("METABRIC", "AGE_AT_DIAGNOSIS", "metabric", "CLAUDIN_SUBTYPE")]:
    clin = pd.read_csv(f"{RAW}/{tag}_clinical.csv", index_col=0)
    clin = clin[~clin.index.duplicated(keep="first")]
    expr = pd.read_csv(f"{RAW}/{tag}_expression_hv.csv", index_col=0)
    panel = pd.read_csv(f"{RAW}/{tag}_expression_sdds_genes.csv", index_col=0)
    tme = pd.read_csv(f"{RAW}/{tag}_expression_tme.csv", index_col=0)
    expr = pd.concat([expr, panel, tme])
    expr = expr[~expr.index.duplicated(keep="first")]
    if tag == "TCGA":
        expr.columns = [c[:12] for c in expr.columns]
        expr = expr.T.groupby(level=0).mean().T
    common = expr.columns.intersection(clin.index)
    expr, clin = expr[common], clin.loc[common]
    os_m = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
    os_s = clin["OS_STATUS"].astype(str).str.startswith("1").astype(int)
    keep = os_m.notna() & (os_m > 0)
    expr, clin = expr.loc[:, keep], clin[keep]
    os_m, os_s = os_m[keep], os_s[keep]

    fr = deconvolute(expr)
    fr.to_csv(f"{OUT}/deconv_fractions_{tag}.csv")

    # subtype-wise z for module score + gene states
    sub = clin[sub_col] if sub_col else None
    if sub is None:
        Zp = expr.loc[["ESR1", "ERBB2"]]
        Zp = Zp.sub(Zp.mean(axis=1), axis=0).div(Zp.std(axis=1) + 1e-9, axis=0)
        sub = pd.Series("TNBC-like", index=expr.columns)
        sub[Zp.loc["ERBB2"] > 1] = "HER2+"
        sub[(Zp.loc["ESR1"] > 0.5) & (sub == "TNBC-like")] = "HR+"
    sub = sub.reindex(expr.columns).fillna("Unknown")
    Z = pd.DataFrame(index=expr.index, columns=expr.columns, dtype=float)
    for s_ in sub.unique():
        cols = sub.index[sub == s_]
        if len(cols) < 5:
            cols = expr.columns
        E = expr[cols]
        Z.loc[:, cols] = E.sub(E.mean(axis=1), axis=0).div(E.std(axis=1) + 1e-9, axis=0)

    esv_mod = Z.loc[[g for g in esv_module if g in Z.index]].mean(axis=0)
    tme_w = ((fr["myCAF"] + fr["SPP1_TAM"]) / 2)
    esv_decon = esv_mod * (tme_w - tme_w.min()) / (tme_w.max() - tme_w.min() + 1e-9)

    # multivariate Cox with deconvolution ESV
    if tag == "TCGA":
        stage_n = clin["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.extract(
            r"Stage ([IV]+)")[0].map({"I": 1, "II": 2, "III": 3, "IV": 4})
    else:
        stage_n = pd.to_numeric(clin["TUMOR_STAGE"], errors="coerce")
    age = pd.to_numeric(clin[age_col], errors="coerce")
    sub_d = pd.get_dummies(sub, prefix="sub", drop_first=True).astype(float)
    mv = pd.DataFrame({"T": os_m, "E": os_s, "ESV_decon": esv_decon,
                       "myCAF_frac": fr["myCAF"], "SPP1TAM_frac": fr["SPP1_TAM"],
                       "age": age, "stage": stage_n})
    mv = pd.concat([mv, sub_d], axis=1).dropna()
    cph = CoxPHFitter()
    cph.fit(mv, "T", "E")
    cph.summary.to_csv(f"{OUT}/cox_multivariate_{tag}_deconv.csv")
    print(f"\n{tag} multivariate (deconv ESV) C-index:",
          round(cph.concordance_index_, 3), flush=True)
    print(cph.summary[["coef", "p"]].round(4).to_string(), flush=True)

    # KM with deconv ESV
    med = esv_decon.median()
    lr = logrank_test(os_m[esv_decon > med], os_m[esv_decon <= med],
                      event_observed_A=os_s[esv_decon > med],
                      event_observed_B=os_s[esv_decon <= med])
    print(tag, "ESV_decon KM log-rank P=", f"{lr.p_value:.2e}", flush=True)
    pd.DataFrame({"T": os_m, "E": os_s,
                  "ESV_decon_high": (esv_decon > med).astype(int)}).assign(
        cohort=tag).to_csv(f"{OUT}/km_data_deconv_{tag}.csv")
print("P0_2_DONE")
