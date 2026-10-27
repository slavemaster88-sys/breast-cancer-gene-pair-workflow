# v3.0.0 — Public reproducibility release

This release accompanies the manuscript **“Reconstruction and exploratory evaluation of rank-based tumour microenvironment-associated gene-pair scores in breast cancer.”** It provides the author-developed analysis code, original manuscript figures, portable provenance manifests, and non-record-level aggregate and quality-assurance results needed for technical inspection of the reconstructed workflow.

## Highlights

- Reconstructed R analysis workflow with archived, privacy-redacted Python provenance scripts.
- Original Figures 1–6 in publication-ready PNG/PDF formats.
- Aggregate drug-model, spatial-section, pathological-response, survival, and audit outputs.
- Portable input locators, source-file hashes, figure-source mappings, software versions, parameters, random seeds, and a self-excluding SHA-256 file manifest.
- Figure 1 provenance documented as author-generated from public GSE176078 count matrices; no published panel was reproduced or adapted.
- Submission-facing terminology corrected and made auditable:
  - `SpaGI-BC` is treated only as a historical internal project/version label.
  - The undefined legacy identifier `sDDS` is described as a **predicted probability of the relative assay-defined sensitive class**.
  - `DOP` is expanded as **durvalumab–olaparib–paclitaxel**.

## Deliberate exclusions

This public release does **not** redistribute frozen model objects, patient-, sample-, cell-, or spot-level derivatives, GDSC-derived response records, substantial KEGG or CellChatDB derivatives, or cBioPortal/METABRIC/TCGA record-level tables. These materials remain in the authors' controlled internal audit archive because they contain source-derived information governed by dataset-specific access or reuse terms. See `THIRD_PARTY_NOTICES.md` and `PROVENANCE_AND_TERMINOLOGY.md`.

## Interpretation

The workflow and outputs support reproducible exploratory evaluation. They do not validate dual-gene mechanisms, treatment-specific benefit, independent clinical performance, or clinical utility. The software and model interfaces are **not validated for clinical use**.

## Licensing and citation

Author-original code is released under the MIT License. Author-original documentation, figures, and eligible aggregate tables are released under CC BY 4.0. Third-party materials are excluded from those grants as specified in `THIRD_PARTY_NOTICES.md`.

Use `CITATION.cff` for software citation metadata. Cite the accompanying manuscript and the archived version 3.0.0 release at https://doi.org/10.5281/zenodo.23684869.
