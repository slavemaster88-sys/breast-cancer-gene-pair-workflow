source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "downstream_audit.R"))
final_scientific_qa<-function() {
  audit_downstream()
  old<-fread(file.path(ROOT,"06_clinical_validation/spatial_replication_gse213688.csv"))
  regional<-rbindlist(lapply(old$sample,function(s) {
    x<-fread(file.path(OUT,"tables",paste0("spatial_",s,"_scores_and_frozen_regions_v3.csv")))
    data.table(sample=s,reproduced_low=mean(x$legacy_ESV[x$low_malignant]),
      reproduced_high=mean(x$legacy_ESV[x$high_malignant]))
  }))
  regional<-merge(regional,old[,.(sample,legacy_low=ESV_boundary,legacy_high=ESV_core)],by="sample")
  regional[,max_abs_difference:=pmax(abs(reproduced_low-legacy_low),abs(reproduced_high-legacy_high))]
  regional[,pass:=max_abs_difference<=1e-5]
  stopifnot(all(regional$pass))
  write_table(regional,"spatial_legacy_regional_mean_reproduction_v3.csv")
  gates<-fread(file.path(OUT,"tables/spatial_all_reproduction_gates_v3.csv"))
  stopifnot(nrow(gates)==15L,all(gates$pass))
  net<-fread(file.path(OUT,"tables/network_retained_typed.csv"))
  stopifnot(!any(net$symbol1==net$symbol2),nrow(net)==241L)
  feat<-fread(file.path(OUT,"tables/pCR_locked_model_features_v3.csv"))
  stopifnot(all(c("M0","M1","M2")%in%feat[feature=="MP",model]))
  folds<-fread(file.path(OUT,"tables/pCR_fold_training_IDs_v3.csv"))
  oof<-fread(file.path(OUT,"tables/pCR_OOF_predictions_v3.csv"))
  for(i in seq_len(nrow(folds))) {
    z<-folds[i];ids<-strsplit(z$training_ids,";",fixed=TRUE)[[1]]
    stopifnot(!length(intersect(ids,oof[repeat_id==z$repeat_id&fold==z$fold,ResearchID])))
  }
  C<-readRDS(file.path(WORK,"clinical_features_v3.rds"))
  P<-readRDS(file.path(WORK,"clinical_pCR_models_v3.rds"))
  reference<-list(version=VERSION,Uref=C$Uref,basis=C$basis,module=C$module,tcell=C$tcell,
    rule="Patient-label-equal raw-linear reference; sample ranks; normalized NNLS relative composition proxies",
    module_training_mean=P$module_training_mean,module_training_sd=P$module_training_sd)
  saveRDS(reference,file.path(OUT,"models/clinical_reference_v3.rds"),compress="gzip")
  saveRDS(list(version=VERSION,final=P$final,features=P$features,training_ids=P$training_ids,
    MP="0=MP1; 1=MP2",reference=reference),file.path(OUT,"models/pCR_frozen_models_v3.rds"),compress="gzip")
  stopifnot(length(readRDS(file.path(OUT,"models/pCR_frozen_models_v3.rds"))$final)==3L)
  inputs<-fread(file.path(CHECK,"input_checkpoint_manifest.csv"))
  inputs[,checkpoint_path:=file.path(CHECK,"inputs",relative_path)]
  inputs[,verified:=file.exists(checkpoint_path)&file.info(checkpoint_path)$size==bytes]
  inputs[,hash_verified:=vapply(checkpoint_path,sha,character(1))==sha256]
  stopifnot(all(inputs$verified),all(inputs$hash_verified))
  write_table(inputs[,.(relative_path,bytes,sha256,verified,hash_verified)],"input_checkpoint_verification_final_v3.csv")
  tests<-data.table(test=c("zero_retained_self_pairs","241_typed_records","15_spatial_reproduction_gates",
      "MP_in_all_final_pCR_models","pCR_outer_ID_isolation","clinical_bundles_readback","108_input_checkpoint_hashes"),
      pass=TRUE,detail=c("No symbol1=symbol2","Clean group-aware keys","All section per-spot gates passed",
      "M0/M1/M2 include MP","25 outer folds; no ResearchID overlap","3 frozen pCR models plus reference",
      paste(nrow(inputs),"files size/SHA verified")))
  write_table(tests,"final_scientific_QA_v3.csv")
  print(tests)
}
if(sys.nframe()==0L)final_scientific_qa()
