#!/usr/bin/env python
"""SpaGI-BC Phase 2, step 4: topological specificity (topS) + network similarity.

For each subtype/TME network (propR-retained edges, weight = rho):
  - node strength w(i) = sum of rho over direct neighbors
  - null: 1,000 permutations of edge weights across the fixed topology
  - topS(i) = (w(i) - mu_R(i)) / sigma_R(i); hubs = topS > 2
Cross-subtype similarity: Spearman of edge weights on shared edges and of
node strength on shared nodes (TNBC vs HER2+ vs ER+).
"""
import itertools
import numpy as np
import pandas as pd
import scipy.stats as st

OUT = "/workspace/SpaGI_BC/02_GI_networks"
N_PERM = 1000
SEED = 0

edges = pd.read_csv(f"{OUT}/propR_edges_all.csv")
edges = edges[edges["retain_propr"]].copy()
ev = pd.read_csv(f"{OUT}/functional_evidence_all.csv", low_memory=False)[
    ["group", "symbol1", "symbol2", "pair_type", "depmap_sl", "depmap_sv",
     "icraft_pass", "functional_pass", "expr_dep_c1", "expr_dep_p1",
     "expr_dep_c2", "expr_dep_p2", "mean_eff_A", "mean_eff_B"]]
ev = ev.drop_duplicates(["group", "symbol1", "symbol2", "pair_type"], keep="first")
edges = edges.merge(ev, on=["group", "symbol1", "symbol2", "pair_type"], how="left")

rng = np.random.default_rng(SEED)
tops_all, hub_rows = [], []
for g, gp in edges.groupby("group"):
    nodes = sorted(set(gp["symbol1"]) | set(gp["symbol2"]))
    nidx = {x: k for k, x in enumerate(nodes)}
    W = np.zeros(len(nodes))
    i = gp["symbol1"].map(nidx).values
    j = gp["symbol2"].map(nidx).values
    w = gp["rho"].values
    np.add.at(W, i, w)
    np.add.at(W, j, w)
    # null: permute edge weights over the fixed topology
    Wnull = np.zeros((N_PERM, len(nodes)))
    for b in range(N_PERM):
        wp = rng.permutation(w)
        np.add.at(Wnull[b], i, wp)
        np.add.at(Wnull[b], j, wp)
    mu, sd = Wnull.mean(0), Wnull.std(0, ddof=1)
    sd[sd == 0] = 1e-12
    topS = (W - mu) / sd
    t = pd.DataFrame({"group": g, "gene": nodes, "strength": W,
                      "topS": topS, "hub": topS > 2})
    t = t.sort_values("topS", ascending=False)
    tops_all.append(t)
    hub_rows.append(t[t["hub"]][["gene", "strength", "topS"]].assign(group=g))
    print(f"{g}: {len(gp)} edges, {len(nodes)} nodes, "
          f"hubs(topS>2) {int(t['hub'].sum())}", flush=True)

tops = pd.concat(tops_all, ignore_index=True)
tops.to_csv(f"{OUT}/topS_all.csv", index=False)
hubs = pd.concat(hub_rows, ignore_index=True)
hubs.to_csv(f"{OUT}/hub_genes.csv", index=False)

# ---- cross-subtype similarity (malignant networks) ----
mal = ["Malignant_TNBC", "Malignant_HER2+", "Malignant_ER+"]
sub = {g: edges[edges["group"] == g] for g in mal}
sim_rows = []
for a, b in itertools.combinations(mal, 2):
    ea, eb = sub[a], sub[b]
    key_a = set(zip(ea["symbol1"], ea["symbol2"], ea["pair_type"]))
    key_b = set(zip(eb["symbol1"], eb["symbol2"], eb["pair_type"]))
    shared = key_a & key_b
    # edge-level Spearman on shared edges
    if len(shared) > 2:
        ma = ea.apply(lambda r: (r["symbol1"], r["symbol2"], r["pair_type"]),
                      axis=1)
        ra = {k: v for k, v in zip(ma, ea["rho"])}
        rb = {k: v for k, v in zip(
            eb.apply(lambda r: (r["symbol1"], r["symbol2"], r["pair_type"]),
                     axis=1), eb["rho"])}
        ks = sorted(shared)
        rho_sp = st.spearmanr([ra[k] for k in ks], [rb[k] for k in ks])
    else:
        rho_sp = (np.nan, np.nan)
    # node-level Spearman on shared nodes
    ta = tops[tops["group"] == a].set_index("gene")["strength"]
    tb = tops[tops["group"] == b].set_index("gene")["strength"]
    common = ta.index.intersection(tb.index)
    node_sp = st.spearmanr(ta[common], tb[common]) if len(common) > 2 else (np.nan, np.nan)
    sim_rows.append({
        "network_a": a, "network_b": b,
        "edges_a": len(ea), "edges_b": len(eb), "shared_edges": len(shared),
        "jaccard_edges": len(shared) / len(key_a | key_b),
        "spearman_edge_rho": rho_sp[0], "spearman_edge_p": rho_sp[1],
        "shared_nodes": len(common),
        "spearman_node_rho": node_sp[0], "spearman_node_p": node_sp[1]})
sim = pd.DataFrame(sim_rows)
sim.to_csv(f"{OUT}/network_similarity.csv", index=False)
print(sim.to_string(index=False))
print("PHASE2_STEP4_DONE", flush=True)
