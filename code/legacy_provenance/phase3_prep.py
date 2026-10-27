#!/usr/bin/env python
"""SpaGI-BC Phase 3, step 1: prepare malignant + TME subsets for CellChat.

TME subsets defined from Phase-1 module scores (obs columns tme_*):
  - FAP+ myCAF: Fibroblasts with high tme_myCAF score (top 50%)
  - SPP1+ TAM:  Myeloid with high tme_SPP1_TAM score (top 50%)
  - CD8+ Tex:   T cells with high tme_CD8_Tex score (top 50%)
  - Endothelial: all Endothelial cells
Malignant epithelial cells by subtype.
Exports: counts (genes x cells) as RDS + metadata CSV, per subtype and pooled.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import pyreadr
import warnings
warnings.filterwarnings("ignore")

ATLAS = "/workspace/SpaGI_BC/01_scRNA_processed/BC_scRNA_annotated.h5ad"
OUT = "/workspace/SpaGI_BC/03_TME_extrinsic_SV"
import os
os.makedirs(OUT, exist_ok=True)

adata = sc.read_h5ad(ATLAS)
obs = adata.obs.copy()

# TME subset assignment
tme_mycaf = obs["tme_myCAF"].values
tme_spp1 = obs["tme_SPP1_TAM"].values
tme_tex = obs["tme_CD8_Tex"].values
ct = obs["celltype"].values

subset = np.full(len(obs), "", dtype=object)
fib = ct == "Fibroblasts"
mye = ct == "Myeloid"
tcell = ct == "T cells"
endo = ct == "Endothelial"
mal = obs["malignancy"].values == "Malignant"

thr_f = np.nanmedian(tme_mycaf[fib])
thr_s = np.nanmedian(tme_spp1[mye])
thr_t = np.nanmedian(tme_tex[tcell])
subset[fib & (tme_mycaf >= thr_f)] = "myCAF"
subset[mye & (tme_spp1 >= thr_s)] = "SPP1_TAM"
subset[tcell & (tme_tex >= thr_t)] = "CD8_Tex"
subset[endo] = "Endothelial"
subset[mal] = "Malignant_" + obs.loc[mal, "subtype"].astype(str).values

keep = subset != ""
print("subset counts:", pd.Series(subset[keep]).value_counts().to_dict(), flush=True)

sub = adata[keep].copy()
sub.obs["cc_subset"] = subset[keep]
# cap each subset at 3000 cells for CellChat speed (random, seed 0)
rng = np.random.default_rng(0)
idx_keep = []
for s, grp in sub.obs.groupby("cc_subset").groups.items():
    gi = np.array(grp)
    if len(gi) > 3000:
        gi = rng.choice(gi, 3000, replace=False)
    idx_keep.extend(gi)
sub = sub[idx_keep].copy()
print("after cap:", sub.obs["cc_subset"].value_counts().to_dict(), flush=True)

X = sp.csr_matrix(sub.raw.X)
counts = pd.DataFrame(np.expm1(X.T).toarray(),
                      index=sub.raw.var_names, columns=sub.obs_names)
# restrict to CellChatDB genes (all CellChat needs) - avoids 4.5GB RDS and
# pyreadr rowname loss
db_genes = set(open("/workspace/cellchatdb_genes.txt").read().split())
counts = counts[counts.index.isin(db_genes)]
# pad zero rows for DB genes absent from the atlas (CellChat requires all
# complex subunits present)
missing = sorted(db_genes - set(counts.index))
pad = pd.DataFrame(0.0, index=missing, columns=counts.columns)
counts = pd.concat([counts, pad])
print("DB genes retained:", counts.shape[0], "padded:", len(missing), flush=True)
counts.reset_index().rename(columns={"index": "gene"}).to_csv(
    f"{OUT}/cellchat_counts.csv.gz", index=False)
meta = sub.obs[["cc_subset", "subtype", "patient", "celltype"]].copy()
meta.to_csv(f"{OUT}/cellchat_meta.csv")
print("exported", counts.shape, flush=True)
print("PHASE3_PREP_DONE", flush=True)
