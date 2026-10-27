#!/usr/bin/env python
"""SpaGI-BC Phase 2, step 3 (v2): functional perturbation evidence for ALL
candidate CGIdb pairs + hard filter for propR-retained edges.

Malignant pairs: DepMap 24Q2 CRISPR gene-effect + CCLE expression in breast
cancer lines (47 lines with both datasets). TusoPerturb dual-perturbation
embeddings are not available in this environment, so we use the standard
single-gene DepMap SL signature:
  - SL evidence: spearman(expr_A, effect_B) or spearman(expr_B, effect_A)
    < -0.2 with P < 0.05 (high expression of one partner renders the other
    essential -> combined loss of viability).
  - SV evidence: positive selection -- spearman > +0.2 (P<0.05) either
    direction, or mean gene effect > +0.05 (KO increases fitness).
DepMap scores: negative = essential (inverted convention already applied:
we correlate expression with the RAW effect score, where more negative =
more essential, so SL = negative correlation).

TME pairs: ICRAFT (Luo et al., Immunity 2025) has no bulk API; we use the
target lists published in the paper: named hits plus their stated pathways
(NF-kB, JAK-STAT, MHC-I antigen presentation, IFN-gamma, TCR signaling,
OXPHOS, DAG kinases).
"""
import numpy as np
import pandas as pd
import scipy.stats as st
import warnings
warnings.filterwarnings("ignore")

DM = "/workspace/SpaGI_BC/00_raw_data/DepMap"
OUT = "/workspace/SpaGI_BC/02_GI_networks"

ICRAFT_CORE = {
    "TNFAIP3", "TRAF2", "BIRC2", "PTPN2", "SOCS1", "TMEM127", "SUSD6",
    "JAK1", "STAT1", "IFNGR1", "IFNGR2", "B2M", "TAP1", "TAP2",
    "CBLB", "CD5", "TMEM222", "ZC3H12A", "CYC1", "DGKA", "DGKZ",
}
ICRAFT_PATHWAYS = {
    # NF-kB, JAK-STAT, MHC-I, IFN-gamma, TCR, OXPHOS, DAG (paper-named pathways)
    "NFKB1", "RELA", "RELB", "NFKB2", "BIRC3", "TRAF6", "IKBKB",
    "JAK2", "STAT2", "STAT3", "STAT5A", "STAT5B", "SOCS3", "PTPN1", "PTPN6",
    "TAPBP", "CALR", "HLA-A", "HLA-B", "HLA-C", "PSMB8", "PSMB9",
    "IFNGR1", "IFNGR2", "IFNG", "CXCL9", "CXCL10", "CCL5",
    "CD3D", "CD3E", "CD3G", "ZAP70", "LCK", "FYN", "CD28", "CTLA4",
    "PDCD1", "CD274", "LAG3", "TIGIT", "HAVCR2", "TOX", "NR4A1",
    "NDUFS1", "NDUFV1", "ATP5F1A", "COX5A", "UQCRC2",
}
ICRAFT_SET = ICRAFT_CORE | ICRAFT_PATHWAYS

print("Loading DepMap...", flush=True)
model = pd.read_csv(f"{DM}/Model.csv")
bres = set(model[model["OncotreeLineage"] == "Breast"]["ModelID"])
effect = pd.read_csv(f"{DM}/CRISPRGeneEffect.csv", index_col=0)
effect = effect.loc[effect.index.isin(bres)]
expr = pd.read_csv(f"{DM}/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
common = effect.index.intersection(expr.index)
effect, expr = effect.loc[common], expr.loc[common]
print(f"breast lines with both: {len(common)}", flush=True)

eff_genes = {c.split(" ")[0]: c for c in effect.columns}
exp_genes = {c.split(" ")[0]: c for c in expr.columns}

pairs = pd.read_csv(f"{OUT}/CGIdb_filtered_pairs.csv").dropna(subset=["symbol1", "symbol2"])
prop = pd.read_csv(f"{OUT}/propR_edges_all.csv")[
    ["symbol1", "symbol2", "pair_type", "rho", "p_perm",
     "pearson_clr", "p_pearson", "retain_propr", "retain_pearson"]]
pairs = pairs.merge(prop, on=["symbol1", "symbol2", "pair_type"], how="left")
pairs["retain_propr"] = pairs["retain_propr"].fillna(False)
print("candidate pairs:", len(pairs), flush=True)

# precompute per-gene vectors for genes appearing in pairs
need = sorted(set(pairs["symbol1"]) | set(pairs["symbol2"]))
eff_vec = {g: effect[eff_genes[g]].values for g in need if g in eff_genes}
exp_vec = {g: expr[exp_genes[g]].values for g in need if g in exp_genes}
print(f"screened genes: {len(eff_vec)}, expressed genes: {len(exp_vec)}", flush=True)

s1 = pairs["symbol1"].values
s2 = pairs["symbol2"].values
n = len(pairs)
c1 = np.full(n, np.nan); p1 = np.full(n, np.nan)
c2 = np.full(n, np.nan); p2 = np.full(n, np.nan)
mA = np.full(n, np.nan); mB = np.full(n, np.nan)
for k in range(n):
    A, B = s1[k], s2[k]
    eA, eB = eff_vec.get(A), eff_vec.get(B)
    if eA is None or eB is None:
        continue
    mA[k], mB[k] = eA.mean(), eB.mean()
    xA, xB = exp_vec.get(A), exp_vec.get(B)
    if xA is not None:
        c1[k], p1[k] = st.spearmanr(xA, eB)
    if xB is not None:
        c2[k], p2[k] = st.spearmanr(xB, eA)

pairs["expr_dep_c1"], pairs["expr_dep_p1"] = c1, p1
pairs["expr_dep_c2"], pairs["expr_dep_p2"] = c2, p2
pairs["mean_eff_A"], pairs["mean_eff_B"] = mA, mB

# functional evidence flags (all candidate pairs)
with np.errstate(invalid="ignore"):
    sl_sig = ((c1 < -0.2) & (p1 < 0.05)) | ((c2 < -0.2) & (p2 < 0.05))
    sv_sig = ((c1 > 0.2) & (p1 < 0.05)) | ((c2 > 0.2) & (p2 < 0.05)) | \
             (mA > 0.05) | (mB > 0.05)
pairs["depmap_sl"] = np.where(sl_sig, True, False)
pairs["depmap_sv"] = np.where(sv_sig, True, False)
pairs["depmap_pass"] = np.where(pairs["pair_type"] == "SL", pairs["depmap_sl"],
                                pairs["depmap_sv"])
in_icraft = pairs["symbol1"].isin(ICRAFT_SET) | pairs["symbol2"].isin(ICRAFT_SET)
pairs["icraft_pass"] = np.where(pairs["group"].str.startswith("Malignant"),
                                False, in_icraft)
pairs["functional_pass"] = np.where(pairs["group"].str.startswith("Malignant"),
                                    pairs["depmap_pass"], pairs["icraft_pass"])

pairs.to_csv(f"{OUT}/functional_evidence_all.csv", index=False)

ret = pairs[pairs["retain_propr"]] if "retain_propr" in pairs.columns else None
print("\nFunctional evidence among ALL candidates:")
print(pairs.groupby(["group", "pair_type"])[["depmap_pass", "icraft_pass"]].sum())
if ret is not None:
    print("\nFunctional pass among propR-retained edges:")
    print(ret.groupby(["group", "pair_type"])["functional_pass"].sum())
print("PHASE2_STEP3_DONE", flush=True)
