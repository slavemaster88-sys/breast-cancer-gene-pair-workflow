# Provenance and terminology decisions

Decision date: 2026-10-27

## Internal project label

Exact searches of PubMed, DOI registries, GitHub, Zenodo, bioRxiv, medRxiv, arXiv, OSF, and Figshare did not identify a formal method or repository corresponding to `SpaGI-BC`. The archive itself uses the term only as a project/version label and provides no expansion. Public-facing manuscript text therefore uses **author-developed breast-cancer gene-pair workflow** instead.

The formally published 2018 SPAGI method (Signaling Pathway Analysis for putative Gene regulatory network Identification; https://doi.org/10.1186/s12918-018-0655-x) is mechanistically and computationally different and is not cited as the source of this workflow.

## Legacy `sDDS` identifier

No formal expansion was found in code, documentation, historical manuscript files, model metadata, or public searches. The prediction interface returns `probability_sensitive`; the documented interpretation is a probability for the relative in-vitro assay-defined sensitive class. Public-facing outputs therefore avoid the undefined acronym and use **predicted sensitive-class probability**. Legacy identifiers remain unchanged only where renaming would compromise reproducibility.

## DOP regimen

The I-SPY2 trial report formally defines DOP as **durvalumab–olaparib–paclitaxel** (Pusztai et al., *Cancer Cell*, 2021; https://doi.org/10.1016/j.ccell.2021.05.009). Figure 6 and the manuscript use the full regimen name. The local code filters the source arm labelled `durvalumab/olaparib`; paclitaxel was the treatment backbone described in the trial report.

## Figure 1

The exact controlled-archive `code/legacy_provenance/phase1_pipeline.py` reads GSE176078 count matrices, performs the archived upstream analysis, assembles a four-panel figure, and writes `Fig1_scRNA_Atlas.pdf/png`. The public-script copy externalises the 26-label subtype mapping so that no record-level mapping is redistributed. The v3 figure script only carries the author-generated figure forward. The figure differs from Wu et al. Figure 1 in panel structure and analysed cell count. It is therefore not a reproduction or adaptation of the published figure.

The authors confirmed on 2026-10-27 that the archived script and figure were created by their team and that the figure was not previously published or licensed to a third party. No third-party figure permission is therefore required; the public source data and original study remain cited.

