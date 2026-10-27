source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
files<-c(
  "01_scRNA_processed/BC_scRNA_annotated.h5ad",
  paste0("02_GI_networks/",c("CGIdb_filtered_pairs.csv","propR_edges_all.csv","functional_evidence_all.csv")),
  paste0("03_TME_extrinsic_SV/",c("cellchat_all_LR.csv","cellchat_counts.csv.gz","cellchat_meta.csv")),
  paste0("00_raw_data/DepMap/",c("Model.csv","OmicsExpressionProteinCodingGenesTPMLogp1.csv")),
  paste0("04_sDDS_model/",c("gdsc_response.csv","prism_response.csv","SCP542_BC_counts.csv",
                            "SCP542_BC_meta.csv","GSE173634_pseudobulk_ensembl.csv")),
  "05_spatial_niche/signature_genes.json",
  "06_clinical_validation/deconv_reference_profiles.csv")
files<-c(files,paste0("00_raw_data/DrugResponse/",
  list.files(file.path(ROOT,"00_raw_data/DrugResponse"),full.names=FALSE,
             pattern="^(GDSC[12]_fitted|PRISM_primary_LFC|PRISM_treatment_info|PRISM_cell_line_info)\\.csv$")))
files<-unique(c(files,
  paste0("00_raw_data/Clinical/",list.files(file.path(ROOT,"00_raw_data/Clinical"))),
  paste0("00_raw_data/Spatial_extra/extracted/",
         list.files(file.path(ROOT,"00_raw_data/Spatial_extra/extracted"),
                    pattern="(filtered_feature_bc_matrix.h5|tissue_positions_list.csv.gz)$")),
  paste0("06_clinical_validation/",list.files(file.path(ROOT,"06_clinical_validation"),
                    pattern="^(spot_ESV_M.*|spatial_replication_gse213688|moransI_gse213688)\\.csv$"))))
sources<-file.path(ROOT,files)
extra<-c(Sys.getenv("SPAGI_GENCODE",file.path(getwd(),"data","reference","gencode.v44.annotation.gtf.gz")),Sys.getenv("SPAGI_KEGG",file.path(getwd(),"data","reference","kegg_2021.json")),
         file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase2_propr.py"),file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase2_step3.py"),file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase2_step4.py"),
         file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase3_tripartite.py"),file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase3_prep.py"),file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase6g_spatial_rep.py"),
         file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase6c_deconv.py"),file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"legacy_provenance","phase1_pipeline.py"))
sources<-c(sources,extra)
rel<-c(files,paste0("original_code_and_reference/",basename(extra)))
stopifnot(all(file.exists(sources)))
m<-data.table(source=sources,relative_path=rel,bytes=file.info(sources)$size)
cat("Input files:",nrow(m),"total GB:",sum(m$bytes)/1e9,"\n")
for(i in seq_len(nrow(m))) {
  dest<-file.path(CHECK,"inputs",m$relative_path[i])
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  stopifnot(system2("cp",c(shQuote(m$source[i]),shQuote(dest)))==0)
  m[i,sha256:=sha(source)]
}
fwrite(m,file.path(CHECK,"input_checkpoint_manifest.csv"))
fwrite(m,file.path(OUT,"input_manifest_v3.csv"))
cat("Copied files; verification is a separate required action.\n")
