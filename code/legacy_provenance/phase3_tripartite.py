#!/usr/bin/env python
"""SpaGI-BC Phase 3, steps 2-3 (memory-safe): tripartite network + ESV score."""
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

OUT = "/workspace/SpaGI_BC/03_TME_extrinsic_SV"
NET = "/workspace/SpaGI_BC/02_GI_networks"
ATLAS = "/workspace/SpaGI_BC/01_scRNA_processed/BC_scRNA_annotated.h5ad"

lr = pd.read_csv(f"{OUT}/cellchat_all_LR.csv")
tme = ["myCAF", "SPP1_TAM", "CD8_Tex", "Endothelial"]
lr = lr[(lr["pval"] < 0.01) & lr["source"].isin(tme) &
        lr["target"].str.startswith("Malignant")].copy()
print("TME->Malignant pairs at P<0.01:", len(lr), flush=True)

edges = pd.read_csv(f"{NET}/propR_edges_all.csv", low_memory=False)
edges = edges[edges["retain_propr"]]
tops = pd.read_csv(f"{NET}/topS_all.csv")
net_nodes, hub = {}, {}
for g in ["Malignant_TNBC", "Malignant_HER2+", "Malignant_ER+"]:
    e = edges[edges["group"] == g]
    net_nodes[g] = set(e["symbol1"]) | set(e["symbol2"])
    tg = tops[tops["group"] == g]
    sv_genes = set(e[e["pair_type"] == "SV"]["symbol1"]) | \
               set(e[e["pair_type"] == "SV"]["symbol2"])
    hub[g] = set(tg[tg["topS"] >= tg["topS"].quantile(0.9)]["gene"]) | \
             (sv_genes & set(tg[tg["topS"] >= tg["topS"].median()]["gene"]))
    print(g, "hubs:", len(hub[g]), flush=True)

kegg = json.load(open("/workspace/kegg_2021.json"))
gene2path = {}
for p, genes in kegg.items():
    for g in genes:
        gene2path.setdefault(g, set()).add(p)

def receptor_upstream_of_hubs(rec, subtype):
    paths = gene2path.get(rec, set())
    hits = set()
    for p in paths:
        hits |= hub[subtype] & set(kegg[p])
    return len(hits) > 0, hits

counts = pd.read_csv(f"{OUT}/cellchat_counts.csv.gz", index_col="gene")
meta = pd.read_csv(f"{OUT}/cellchat_meta.csv", index_col=0)
markers = {}
for s in ["myCAF", "SPP1_TAM"]:
    own = meta.index[meta["cc_subset"] == s]
    rest = meta.index[meta["cc_subset"] != s]
    mu_own = counts.loc[:, own].mean(axis=1)
    frac_own = (counts.loc[:, own] > 0).mean(axis=1)
    mu_rest = counts.loc[:, rest].mean(axis=1)
    score = mu_own - mu_rest
    markers[s] = set(score[(frac_own > 0.1) & (score > 0)].nlargest(100).index)
    print(s, "markers:", len(markers[s]), flush=True)

rows = []
for _, r in lr.iterrows():
    src, tgt = r["source"], r["target"]
    lig_ok = (r["ligand"] in markers.get(src, set())) or \
             any(r["ligand"] in net_nodes[g] for g in net_nodes)
    rec_node = r["receptor"] in net_nodes.get(tgt, set())
    up, hits = receptor_upstream_of_hubs(r["receptor"], tgt)
    if lig_ok and (rec_node or up):
        rows.append({
            "ligand": r["ligand"], "receptor": r["receptor"],
            "tme_source": src, "malignant_target": tgt,
            "prob": r["prob"], "pval": r["pval"],
            "pathway_name": r["pathway_name"],
            "ligand_is_tme_marker": r["ligand"] in markers.get(src, set()),
            "ligand_is_SL_SV_node": any(r["ligand"] in net_nodes[g]
                                        for g in net_nodes),
            "receptor_is_SL_SV_node": rec_node,
            "receptor_upstream_SV_hubs": ";".join(sorted(hits)) if up else "",
            "n_upstream_SV_hubs": len(hits)})
tri = pd.DataFrame(rows)
tri.to_csv(f"{OUT}/Tripartite_LR_SV_network.csv", index=False)
print("tripartite pairs:", len(tri), flush=True)
print(tri.groupby(["tme_source", "malignant_target"]).size(), flush=True)

# ---------------- ESV score (memory-safe) ----------------
adata = sc.read_h5ad(ATLAS)
obs = adata.obs
X = sp.csr_matrix(adata.raw.X)
var_names = adata.raw.var_names

# all genes we need: module genes + ligands
need = set()
for _, r in tri.iterrows():
    need.add(r["receptor"])
    need.update(str(r["receptor_upstream_SV_hubs"]).split(";"))
    need.add(r["ligand"])
need = {g for g in need if g in set(var_names)}
sel_rows = np.where(var_names.isin(need))[0]
sub = X[:, sel_rows]  # sparse cells x genes
sub.data = np.expm1(sub.data)
E = pd.DataFrame(sub.toarray().T, index=var_names[sel_rows], columns=obs.index)
print("needed-gene matrix:", E.shape, flush=True)

# TME weight per patient: myCAF/SPP1_TAM abundance x mean ligand expression
fib = (obs["celltype"] == "Fibroblasts").values
mye = (obs["celltype"] == "Myeloid").values
sel_tme = np.zeros(len(obs), bool)
sel_tme[fib & (obs["tme_myCAF"].values >= np.nanmedian(obs.loc[fib, "tme_myCAF"]))] = True
sel_tme[mye & (obs["tme_SPP1_TAM"].values >= np.nanmedian(obs.loc[mye, "tme_SPP1_TAM"]))] = True
ligands = sorted(set(tri["ligand"]) & set(E.index))
lig_expr = E.loc[ligands]  # genes x all cells
lig_tme = lig_expr.loc[:, lig_expr.columns[sel_tme]]
lig_pat = lig_tme.T.mean(axis=1).groupby(obs.loc[lig_tme.columns, "patient"].values).mean()
abund = pd.Series(sel_tme).groupby(obs["patient"].values).mean()
weight_pat = lig_pat * abund
print("patient weights:", weight_pat.describe().round(3).to_dict(), flush=True)

mal = (obs["malignancy"] == "Malignant").values
esv_rows = []
for sub_t in ["TNBC", "HER2+", "ER+"]:
    g = f"Malignant_{sub_t}"
    tris = tri[tri["malignant_target"] == g]
    recs = set(tris["receptor"])
    hubs = set()
    for h in tris["receptor_upstream_SV_hubs"].dropna():
        hubs |= set(str(h).split(";"))
    module = sorted((recs | hubs) & set(E.index))
    if not module:
        continue
    M = E.loc[module].T
    Z = (M - M.mean()) / (M.std() + 1e-9)
    module_score = Z.mean(axis=1)
    w = weight_pat.reindex(obs.loc[module_score.index, "patient"].values).fillna(0).values
    out = pd.DataFrame({
        "subtype": obs.loc[module_score.index, "subtype"].astype(str).values,
        "patient": obs.loc[module_score.index, "patient"].values,
        f"module_{sub_t}": module_score.values,
        f"ESV_{sub_t}": module_score.values * w})
    esv_rows.append(out.set_index(module_score.index))
    print(f"ESV_{sub_t}: module {len(module)} genes, ligands {len(ligands)}",
          flush=True)

esv = pd.concat(esv_rows, axis=1)
esv.to_csv(f"{OUT}/ESV_scores_per_cell.csv")
print("PHASE3_TRIPARTITE_DONE", flush=True)
