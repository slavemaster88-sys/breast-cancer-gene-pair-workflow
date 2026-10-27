#!/usr/bin/env python
"""SpaGI-BC Phase 2, step 2: pair-restricted propR (rho) + Pearson benchmark.

rho definition verified from propr source (src/lr2propr.cpp, lr2rho):
    rho_ij = 1 - VLR_ij / (var_i + var_j)
           = 2 * cov(clr_i, clr_j) / (var_i + var_j)
with clr = per-sample centered log counts (ivar="clr").
Empirical p-values: sample-wise permutation of cells, B=200.
"""
import os
import numpy as np
import pandas as pd
import pyreadr
import scipy.stats as st
import warnings
warnings.filterwarnings("ignore")

EXPORT = "/workspace/SpaGI_BC/02_GI_networks/_export"
OUT = "/workspace/SpaGI_BC/02_GI_networks"
B_PERM = 200
RHO_MIN = 0.3
ALPHA = 0.05

GROUPS = ["Malignant_TNBC", "Malignant_HER2+", "Malignant_ER+",
          "T_cells", "Myeloid", "Fibroblasts", "Endothelial"]


def clr_transform(counts):
    """counts: genes x cells. propr-exact: zeros -> min nonzero, log, center
    per sample. Returns (clr, keep_mask) with zero-variance genes dropped."""
    nz = counts[counts > 0]
    ct = counts.copy()
    ct[ct == 0] = nz.min()
    lr = np.log(ct)
    clr = lr - lr.mean(axis=0, keepdims=True)  # center per cell (sample)
    keep = clr.var(axis=1, ddof=1) > 0
    return clr[keep], keep


def rho_pairs(clr, i_idx, j_idx):
    """rho for gene pairs; clr is genes x cells."""
    X, Y = clr[i_idx], clr[j_idx]
    n = clr.shape[1]
    mx = X.mean(axis=1, keepdims=True)
    my = Y.mean(axis=1, keepdims=True)
    cov = ((X - mx) * (Y - my)).sum(axis=1) / (n - 1)
    var_x = X.var(axis=1, ddof=1)
    var_y = Y.var(axis=1, ddof=1)
    return 2 * cov / (var_x + var_y)


def perm_pvalues(clr, i_idx, j_idx, obs_rho, B=B_PERM, seed=0):
    """Empirical one-sided p: P(rho_perm >= rho_obs).

    Null: no association between gene i and gene j profiles across cells.
    Each gene's expression profile is shuffled independently across cells
    (per-row permutation), breaking the i-j pairing while preserving each
    gene's marginal distribution. (Shuffling whole cells would preserve the
    correlation and give a degenerate null.)
    """
    rng = np.random.default_rng(seed)
    n_cells = clr.shape[1]
    X, Y = clr[i_idx], clr[j_idx]
    ge = np.zeros(len(obs_rho))
    for b in range(B):
        ix = np.argsort(rng.random(X.shape), axis=1)
        jx = np.argsort(rng.random(Y.shape), axis=1)
        rp = rho_pairs(np.vstack([np.take_along_axis(X, ix, axis=1),
                                  np.take_along_axis(Y, jx, axis=1)]),
                       np.arange(len(obs_rho)),
                       len(obs_rho) + np.arange(len(obs_rho)))
        ge += (rp >= obs_rho)
    return (1 + ge) / (1 + B)


def main():
    pairs = pd.read_csv(f"{OUT}/CGIdb_filtered_pairs.csv")
    pairs = pairs.dropna(subset=["symbol1", "symbol2"])
    out_rows = []
    for g in GROUPS:
        f = f"{EXPORT}/expr_{g}.rds"
        if not os.path.exists(f):
            print("missing", f, flush=True)
            continue
        r = pd.read_parquet(f"{EXPORT}/expr_{g}.parquet").set_index("gene")
        expr = r
        genes = list(expr.index)
        gidx = {x: k for k, x in enumerate(genes)}
        counts = expr.values.astype(np.float64)  # genes x cells
        clr, keepv = clr_transform(counts)
        kept_genes = [genes[k] for k in np.where(keepv)[0]]
        gidx = {x: k for k, x in enumerate(kept_genes)}

        gp = pairs[pairs["group"] == g]
        keep = gp["symbol1"].isin(gidx) & gp["symbol2"].isin(gidx)
        gp = gp[keep]
        i_idx = np.array([gidx[x] for x in gp["symbol1"]])
        j_idx = np.array([gidx[x] for x in gp["symbol2"]])
        print(f"{g}: {len(gp)} candidate pairs, {clr.shape[0]} genes, "
              f"{clr.shape[1]} cells", flush=True)

        rho = rho_pairs(clr, i_idx, j_idx)
        # Pearson baseline on the same clr data (vectorized)
        X, Y = clr[i_idx], clr[j_idx]
        Xc = X - X.mean(axis=1, keepdims=True)
        Yc = Y - Y.mean(axis=1, keepdims=True)
        pear = (Xc * Yc).sum(axis=1) / np.sqrt(
            (Xc ** 2).sum(axis=1) * (Yc ** 2).sum(axis=1))
        dfree = X.shape[1] - 2
        tstat = pear * np.sqrt(dfree / (1 - pear ** 2 + 1e-12))
        ppear = 2 * st.t.sf(np.abs(tstat), dfree)

        # permutation p-values only for pairs passing rho threshold
        cand = rho > RHO_MIN
        pval = np.ones(len(gp))
        if cand.sum() > 0:
            pval[cand] = perm_pvalues(clr, i_idx[cand], j_idx[cand], rho[cand])

        res = gp[["symbol1", "symbol2", "pair_type", "IntegrateScore",
                  "source", "tissue"]].copy()
        res.insert(0, "group", g)
        res["rho"] = rho
        res["p_perm"] = pval
        res["pearson_clr"] = pear
        res["p_pearson"] = ppear
        res["retain_propr"] = (res["rho"] > RHO_MIN) & (res["p_perm"] < ALPHA)
        res["retain_pearson"] = (res["pearson_clr"] > RHO_MIN) & (res["p_pearson"] < ALPHA)
        out_rows.append(res)
        res.to_csv(f"{OUT}/propR_edges_{g}.csv", index=False)
        print(f"  retained: propR {res['retain_propr'].sum()} vs "
              f"Pearson {res['retain_pearson'].sum()} pairs", flush=True)

    allres = pd.concat(out_rows, ignore_index=True)
    allres.to_csv(f"{OUT}/propR_edges_all.csv", index=False)
    print("PHASE2_PROPR_DONE", flush=True)


if __name__ == "__main__":
    main()
