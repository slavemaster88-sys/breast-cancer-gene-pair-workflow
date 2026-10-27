source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "clinical_models.R"))
audit_downstream<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  C<-readRDS(file.path(WORK,"clinical_features_v3.rds"))
  # Temporal coverage gate: the panel and threshold are not relaxed for this source.
  f<-file.path(ROOT,"00_raw_data/GSE114461/GSE114461_Bortezomib.csv.gz")
  g<-fread(f,skip=1,header=FALSE,select=1)[[1]];sy<-map_symbols(g,D$mapping)
  cv<-length(intersect(D$U,sy))/length(D$U)
  gate<-data.table(dataset="GSE114461",n_source_genes=length(g),n_measurable_U=length(intersect(D$U,sy)),
     U_size=length(D$U),coverage=cv,required_coverage=.95,status=if(cv<.95)"rejected_no_v3_prediction" else "eligible",
     consequence="Old temporal score and significance removed; no missing-score zero fill")
  write_table(gate,"temporal_coverage_gate_v3.csv")
  if(cv>=.95)stop("Temporal source passes coverage; explicit scoring implementation required")
  net<-readRDS(file.path(WORK,"network_v3.rds"))
  typed<-net$retained[(symbol1=="CD36"&symbol2=="PCYOX1L")|(symbol1=="PCYOX1L"&symbol2=="CD36")]
  write_table(data.table(pair="CD36__PCYOX1L",n_clean_retained_group_records=nrow(typed),
    v3_status=if(nrow(typed))"clean_network_member" else "historical_only_not_clean_network",
    manuscript_action="Remove v3 selected-pair efficacy claim; retain provenance note without reselection"),
    "historical_CD36_PCYOX1L_status_v3.csv")
  # Composition proxy and clinical transformations must remain single-sample/batch invariant.
  X<-read_microarray("Pusztai2021_geneLevel.txt.gz",D$mapping)
  full<-rank_ref(X[1:4,,drop=FALSE],C$Uref);single<-rank_ref(X[1,,drop=FALSE],C$Uref)
  Pfull<-composition_proxy(full,C$basis);Pone<-composition_proxy(single,C$basis)
  stopifnot(max(abs(Pfull[1,]-Pone[1,]))<1e-12)
  check<-data.table(test="clinical_reference_single_sample_batch_invariance",
    max_abs_difference=max(abs(Pfull[1,]-Pone[1,])),pass=TRUE)
  write_table(check,"clinical_reference_invariance_v3.csv")
  # Paired conditional confidence intervals use aligned ResearchIDs, not cell-level replicates.
  oof<-fread(file.path(OUT,"tables/pCR_OOF_predictions_v3.csv"))
  stopifnot(!anyDuplicated(oof[,.(ResearchID,model,repeat_id)]))
  p<-oof[,.(y=y[1],p=mean(probability_pCR)),by=.(model,ResearchID)]
  pairs<-list(c("M1","M0"),c("M2","M1"),c("M2","M0"))
  diff<-rbindlist(lapply(pairs,function(ms) {
    a<-p[model==ms[1],.(ResearchID,y,a=p)]
    b<-p[model==ms[2],.(ResearchID,b=p)]
    z<-merge(a,b,by="ResearchID");n<-nrow(z);set.seed(SEED)
    bs<-replicate(2000L,{i<-sample.int(n,n,replace=TRUE);auc_rank(z$y[i],z$a[i])-auc_rank(z$y[i],z$b[i])})
    data.table(comparison=paste(ms,collapse=" minus "),n=n,
      AUROC_meanOOF_difference=auc_rank(z$y,z$a)-auc_rank(z$y,z$b),
      conditional_CI_lo=quantile(bs,.025,names=FALSE),conditional_CI_hi=quantile(bs,.975,names=FALSE))
  }))
  write_table(diff,"pCR_paired_conditional_comparisons_v3.csv")
  miss<-rbindlist(lapply(names(C$features),function(tag) {
    f<-C$features[[tag]]
    data.table(cohort=tag,feature=setdiff(names(f),"sample_id"),
      n=nrow(f),n_missing=vapply(f[,-1],function(x)sum(!is.finite(x)),integer(1)))
  }))
  write_table(miss,"clinical_feature_missingness_v3.csv")
  for(n in c("network_v3.rds","modules_v3.rds","training_data_v3.rds","pilot_model_v3.rds","pilot_spatial_v3.rds",
             "cross_assay_expression_v3.rds","atlas_sDDS_v3.rds","spatial_v3.rds","clinical_features_v3.rds",
             "clinical_pCR_models_v3.rds","survival_v3.rds")) {
    a<-file.path(WORK,n);b<-file.path(CHECK,n)
    stopifnot(file.info(b)$size>0,sha(a)==sha(b));invisible(readRDS(b))
  }
  # Additional inputs used for legacy reproduction, temporal gating and new public records.
  old<-fread(file.path(OUT,"input_manifest_v3.csv"))
  extra<-c(file.path(ROOT,"03_TME_extrinsic_SV/Tripartite_LR_SV_network.csv"),f,
           Sys.getenv("SPAGI_SPATIAL_SOFT",file.path(getwd(),"data","reference","GSE213688_family.soft.gz")),
           file.path(CHECK,c("METABRIC_new_module_gene_records_v3.rds","TCGA_new_module_gene_records_v3.rds",
                            "brca_metabric_sample_manifest_v3.rds","brca_tcga_pan_can_atlas_2018_sample_manifest_v3.rds")))
  extra<-extra[file.exists(extra)]
  added<-data.table(source=extra,relative_path=paste0("v3_additions/",basename(extra)),
                   bytes=as.numeric(file.info(extra)$size),sha256=vapply(extra,sha,character(1)))
  for(i in seq_len(nrow(added))) {
    dest<-file.path(CHECK,"inputs",added$relative_path[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    stopifnot(system2("cp",c(shQuote(added$source[i]),shQuote(dest)))==0)
  }
  merged<-rbindlist(list(old[!source%in%added$source],added),use.names=TRUE)
  fwrite(merged,file.path(OUT,"input_manifest_v3.csv"))
  fwrite(merged,file.path(CHECK,"input_checkpoint_manifest.csv"))
  print(gate);print(diff)
  invisible(list(gate=gate,comparisons=diff))
}
if(sys.nframe()==0L)audit_downstream()
