#!/usr/bin/env python
"""SpaGI-BC Phase 1: BC scRNA-seq atlas curation, QC, annotation, CNV-based
malignant cell identification (GSE176078, Wu et al. 2021).
Scanpy + harmonypy + infercnvpy implementation of the DISCERN-style protocol.
"""
import os, glob, tarfile, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import harmonypy

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
sc.settings.figdir = "/workspace/SpaGI_BC/figures"
matplotlib.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]

RAW = "/workspace/SpaGI_BC/00_raw_data/GSE176078"
PROC = "/workspace/SpaGI_BC/01_scRNA_processed"
FIG = "/workspace/SpaGI_BC/figures"
GTF = "/workspace/ref/gencode.v44.annotation.gtf.gz"
CKPT = "/workspace/SpaGI_BC/01_scRNA_processed/_checkpoint_atlas.h5ad"

SUBTYPE3 = lambda s: "HER2+" if "HER2" in s else s

def load_subtype_map():
    """Load public sample labels and clinical subtypes supplied by the user.

    The exact 26-label mapping used for the archived Figure 1 is retained in
    the authors' controlled audit archive and omitted from this public-code
    candidate under its zero-record-level release policy. Prepare a two-column
    CSV (sample_label, subtype) from the authorised GSE176078 source metadata
    and point SPAGI_SUBTYPE_MAP to it.
    """
    path = os.environ.get("SPAGI_SUBTYPE_MAP")
    if not path:
        raise RuntimeError(
            "Set SPAGI_SUBTYPE_MAP to an authorised CSV with columns "
            "sample_label and subtype."
        )
    mapping = pd.read_csv(path)
    required = {"sample_label", "subtype"}
    missing = required.difference(mapping.columns)
    if missing:
        raise ValueError(f"Subtype mapping is missing columns: {sorted(missing)}")
    if mapping["sample_label"].duplicated().any():
        raise ValueError("Subtype mapping contains duplicate sample_label values")
    return dict(zip(mapping["sample_label"], mapping["subtype"]))

# ---------------------------------------------------------------- load
def load_all():
    subtype_map = load_subtype_map()
    adatas = []
    for tar in sorted(glob.glob(f"{RAW}/*.tar.gz")):
        cid = os.path.basename(tar).split("_")[1].replace(".tar.gz", "")
        base = f"{RAW}/extracted/{cid}"
        if not os.path.isdir(base):
            with tarfile.open(tar) as tf:
                tf.extractall(base)
        mtx = glob.glob(f"{base}/**/count_matrix_sparse.mtx*", recursive=True)
        d = os.path.dirname(mtx[0])
        a = sc.read_mtx(f"{d}/count_matrix_sparse.mtx").T  # cells x genes
        genes = pd.read_csv(glob.glob(f"{d}/count_matrix_genes.tsv*",
                                      recursive=True)[0], sep="\t", header=None)
        a.var_names = genes[1].values if genes.shape[1] > 1 else genes[0].values
        a.var["gene_ids"] = genes[0].values
        a.obs_names = pd.read_csv(glob.glob(f"{d}/count_matrix_barcodes.tsv*",
                                            recursive=True)[0],
                                  sep="\t", header=None)[0].values
        # optional per-cell metadata shipped with the sample
        meta_f = glob.glob(f"{base}/**/metadata.csv", recursive=True)
        if meta_f:
            md = pd.read_csv(meta_f[0], index_col=0)
            shared = a.obs_names.intersection(md.index)
            if len(shared):
                for c in md.columns:
                    a.obs.loc[shared, f"w_{c}"] = md.loc[shared, c].values
        a.obs["patient"] = cid
        if cid not in subtype_map:
            raise KeyError(f"No authorised subtype mapping supplied for {cid}")
        a.obs["subtype_clinical"] = subtype_map[cid]
        a.obs["subtype"] = SUBTYPE3(subtype_map[cid])
        a.var_names_make_unique()
        adatas.append(a)
        print(f"{cid}: {a.shape}", flush=True)
    adata = sc.concat(adatas, join="inner")
    adata.obs_names_make_unique()
    print("Combined:", adata.shape, flush=True)
    return adata

# ---------------------------------------------------------------- QC
def run_harmony(adata, basis="X_pca", adjusted_basis="X_harmony", key="patient"):
    """Direct harmonypy call, robust to scanpy/harmonypy API drift."""
    ho = harmonypy.run_harmony(adata.obsm[basis], adata.obs, key, verbose=False)
    Z = np.asarray(ho.Z_corr)
    adata.obsm[adjusted_basis] = Z.T if Z.shape[0] == adata.obsm[basis].shape[1] else Z
    return adata


def qc_and_normalize(adata):
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)          # (a)
    log_counts = np.log10(adata.obs["total_counts"])  # (b) 3-MAD on log10 UMI
    med = log_counts.median()
    mad = (log_counts - med).abs().median()
    keep = (log_counts > med - 3 * mad) & (log_counts < med + 3 * mad)
    adata = adata[keep.values, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= 20, :].copy()  # (c)
    print(f"QC: {n0} -> {adata.n_obs} cells ({n0 - adata.n_obs} removed)", flush=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="patient")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    run_harmony(adata, basis="X_pca", adjusted_basis="X_harmony", key="patient")
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph",
                 n_iterations=2)
    return adata

# ---------------------------------------------------------------- annotation
MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T cells": ["CD3D", "CD3E", "CD2", "TRAC"],
    "B cells": ["CD79A", "MS4A1", "CD79B"],
    "Myeloid": ["CD68", "LYZ", "CD14", "FCGR3A"],
    "Fibroblasts": ["COL1A1", "DCN", "THY1", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Plasma": ["MZB1", "JCHAIN", "XBP1"],
}
TME_MARKERS = {
    "myCAF": ["FAP", "ACTA2", "TAGLN"],
    "iCAF": ["CXCL12", "IL6", "CFD"],
    "SPP1_TAM": ["SPP1", "APOE"],
    "C1QC_TAM": ["C1QC", "C1QB"],
    "CD8_Tex": ["PDCD1", "LAG3", "TIGIT", "TOX"],
    "CD8_Teff": ["GZMB", "PRF1", "NKG7"],
}

def annotate(adata):
    for ct, genes in MARKERS.items():
        gl = [g for g in genes if g in adata.raw.var_names]
        sc.tl.score_genes(adata, gene_list=gl, score_name=f"s_{ct.replace(' ', '_')}",
                    use_raw=True)
    score_cols = [f"s_{ct.replace(' ', '_')}" for ct in MARKERS]
    S = adata.obs[score_cols].values
    adata.obs["celltype"] = [list(MARKERS)[i] for i in S.argmax(axis=1)]
    for pop, genes in TME_MARKERS.items():
        gl = [g for g in genes if g in adata.raw.var_names]
        sc.tl.score_genes(adata, gene_list=gl, score_name=f"tme_{pop}", use_raw=True)
    print(adata.obs["celltype"].value_counts(), flush=True)
    return adata

# ---------------------------------------------------------------- CNV / malignancy
def add_gene_positions(adata, gtf_path):
    """Parse GENCODE GTF and attach chromosome/start/end to adata.var."""
    import gzip
    cols = ["gene_id", "gene_name", "chromosome", "start", "end"]
    rows = []
    op = gzip.open if gtf_path.endswith(".gz") else open
    with op(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            attr = dict(kv.strip().split('"')[0:2][0:1] + [None] if False else [])
            d = {}
            for kv in p[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                k, _, v = kv.partition(" ")
                d[k] = v.strip('"')
            rows.append((d.get("gene_id"), d.get("gene_name"), p[0], int(p[3]), int(p[4])))
    g = pd.DataFrame(rows, columns=cols).drop_duplicates("gene_name")
    g = g.set_index("gene_name")
    common = adata.var_names.intersection(g.index)
    adata.var["chromosome"] = np.nan
    adata.var["start"] = np.nan
    adata.var["end"] = np.nan
    adata.var.loc[common, "chromosome"] = g.loc[common, "chromosome"].values
    adata.var.loc[common, "start"] = g.loc[common, "start"].values
    adata.var.loc[common, "end"] = g.loc[common, "end"].values
    # infercnvpy expects these dtypes
    adata.var["chromosome"] = adata.var["chromosome"].astype("category")
    adata.var["start"] = adata.var["start"].fillna(-1).astype("int64")
    adata.var["end"] = adata.var["end"].fillna(-1).astype("int64")
    print(f"Gene positions added for {len(common)}/{adata.n_vars} genes", flush=True)


def malignant_calling(adata):
    # epithelial re-clustering into 7 groups
    epi = adata[adata.obs["celltype"] == "Epithelial"].copy()
    sc.pp.highly_variable_genes(epi, n_top_genes=3000)
    sc.pp.scale(epi, max_value=10)
    sc.tl.pca(epi, n_comps=30)
    run_harmony(epi, basis="X_pca", adjusted_basis="X_harmony_epi", key="patient")
    sc.pp.neighbors(epi, use_rep="X_harmony_epi", n_neighbors=30)
    sc.tl.umap(epi)
    # scan resolutions, keep the one yielding exactly 7 clusters (or closest)
    trials = []
    for res in np.arange(0.20, 0.61, 0.01):
        sc.tl.leiden(epi, resolution=float(res), key_added="epi_leiden",
                     flavor="igraph", n_iterations=2)
        k = epi.obs["epi_leiden"].nunique()
        trials.append((float(res), k))
        if k == 7:
            break
    res, k = min(trials, key=lambda t: abs(t[1] - 7))
    if k != 7:
        sc.tl.leiden(epi, resolution=res, key_added="epi_leiden",
                     flavor="igraph", n_iterations=2)
        print(f"WARNING: exact 7 clusters not reached; using k={k} at res={res}",
              flush=True)
    print(f"Epithelial clusters (res={res:.2f}):",
          dict(epi.obs["epi_leiden"].value_counts()), flush=True)

    # inferCNV: epithelial + subsampled T cells (diploid reference) together;
    # n_jobs=1 keeps joblib from duplicating the AnnData in worker processes
    rng = np.random.default_rng(0)
    tcells = adata[adata.obs["celltype"] == "T cells"].obs_names
    tref = rng.choice(tcells.values, size=min(3000, len(tcells)), replace=False)
    cnv_adata = adata[
        adata.obs["celltype"].isin(["Epithelial", "T cells"])
        & adata.obs_names.isin(tref) | (adata.obs["celltype"] == "Epithelial")
    ].copy()
    add_gene_positions(cnv_adata, GTF)
    print("CNV input:", cnv_adata.shape, flush=True)
    cnv.tl.infercnv(cnv_adata, reference_key="celltype",
                    reference_cat=["T cells"], window_size=250,
                    n_jobs=1, chunksize=2000)
    cnv.tl.pca(cnv_adata)
    epi_mask = (cnv_adata.obs["celltype"] == "Epithelial").values
    epi.obs["cnv_score"] = np.abs(cnv_adata.obsm["X_cnv"][epi_mask]).mean(axis=1)
    epi.obs["cnv_pca1"] = cnv_adata.obsm["X_cnv_pca"][epi_mask, 0]

    # classify: a cluster is normal if it lies within the T-cell diploid cloud
    # (median distance to T-cell centroid <= 95th pct of T cells' own radii)
    tc_mask = (cnv_adata.obs["celltype"] == "T cells").values
    P = cnv_adata.obsm["X_cnv_pca"]
    centroid = P[tc_mask].mean(axis=0)
    d_all = np.sqrt(((P - centroid) ** 2).sum(axis=1))
    thr = np.percentile(d_all[tc_mask], 95)
    cnv_adata.obs["dist_to_diploid"] = d_all
    epi.obs["dist_to_diploid"] = d_all[~tc_mask]
    cl = epi.obs.groupby("epi_leiden", observed=True)["dist_to_diploid"].median()
    # cell-level rule: cells inside the T-cell diploid cloud are normal epithelial
    epi.obs["malignancy"] = np.where(epi.obs["dist_to_diploid"] <= thr,
                                     "Normal epithelial", "Malignant")
    tmed = epi.obs.groupby("epi_leiden", observed=True)["cnv_score"].median()
    print("Distance-to-diploid medians per cluster:\n", cl.round(3), flush=True)
    print("Diploid radius threshold (95th pct T cells):", round(thr, 3), flush=True)
    print("CNV medians per cluster:\n", tmed, flush=True)
    print("Malignant cells:", int((epi.obs["malignancy"]=="Malignant").sum()),
          "Normal epithelial:", int((epi.obs["malignancy"]=="Normal epithelial").sum()), flush=True)
    return adata, epi

# ---------------------------------------------------------------- figure
def figure1(adata, epi):
    tmp = adata.copy()
    tmp.obs["malignancy"] = "Non-epithelial"
    tmp.obs.loc[epi.obs_names, "malignancy"] = epi.obs["malignancy"].values

    fig, axes = plt.subplots(2, 2, figsize=(17, 13))
    sc.pl.umap(adata, color="celltype", ax=axes[0, 0], show=False,
               frameon=False, title="A  Major cell types")
    sc.pl.umap(adata, color="subtype", ax=axes[0, 1], show=False, frameon=False,
               title="B  Molecular subtype",
               palette={"ER+": "#0279EE", "HER2+": "#FF9400", "TNBC": "#E9ED4C"})
    sc.pl.umap(tmp, color="malignancy", ax=axes[1, 0], show=False, frameon=False,
               title="C  inferCNV malignancy status",
               palette={"Malignant": "#C0392B", "Normal epithelial": "#75A025",
                        "Non-epithelial": "#DDDDDD"})
    ax = axes[1, 1]
    prop = (tmp.obs.groupby(["patient", "celltype"], observed=True)
            .size().unstack(fill_value=0))
    prop = prop.div(prop.sum(axis=1), axis=0)
    order = sorted(prop.index, key=lambda p: (SUBTYPE3(SUBTYPE[p]), p))
    prop.loc[order].plot(kind="bar", stacked=True, ax=ax, width=0.85,
                         colormap="tab20", legend=False)
    ax.set_title("D  Cell type proportions per patient")
    ax.set_xlabel(""); ax.set_ylabel("Fraction")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    for j, p in enumerate(order):
        c = {"ER+": "#0279EE", "HER2+": "#FF9400", "TNBC": "#E9ED4C"}[SUBTYPE3(SUBTYPE[p])]
        ax.axvspan(j - 0.5, j + 0.5, ymin=-0.06, ymax=0, color=c, clip_on=False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=6, loc="center left",
              bbox_to_anchor=(1.0, 0.5), frameon=False)
    plt.tight_layout()
    fig.savefig(f"{FIG}/Fig1_scRNA_Atlas.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{FIG}/Fig1_scRNA_Atlas.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    if os.path.exists(CKPT):
        print("Loading checkpoint...", flush=True)
        adata = sc.read_h5ad(CKPT)
    else:
        adata = load_all()
        adata = qc_and_normalize(adata)
        adata = annotate(adata)
        adata.write(CKPT, compression="gzip")
    adata, epi = malignant_calling(adata)
    for col in ["epi_leiden", "malignancy", "cnv_score", "cnv_pca1"]:
        adata.obs[col] = epi.obs[col]
    adata.uns["phase1"] = {
        "reference": "T cells (diploid)",
        "n_epithelial_clusters": 7,
        "qc": "genes>=200; UMI 3-MAD (log10); mito<=20%",
    }
    adata.write(f"{PROC}/BC_scRNA_annotated.h5ad", compression="gzip")
    epi.write(f"{PROC}/BC_epithelial_cnv.h5ad", compression="gzip")
    figure1(adata, epi)
    summ = (adata.obs.groupby(["celltype", "malignancy"], observed=True)
            .size().rename("n_cells").reset_index())
    summ.to_csv("/workspace/SpaGI_BC/tables/Fig1_cell_counts.csv", index=False)
    print("PHASE1_DONE", flush=True)
