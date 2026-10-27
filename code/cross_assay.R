source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "prepare_training.R"))
source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "train_models.R"))
prepare_cross_assay<-function(D) {
  models<-fread(file.path(ROOT,"00_raw_data/DepMap/Model.csv"))
  normalize<-function(x)gsub("[^A-Z0-9]","",toupper(sub("_BREAST$","",x)))
  lookup<-unique(rbindlist(list(models[,.(name=normalize(CellLineName),model_id=ModelID)],
             models[,.(name=normalize(StrippedCellLineName),model_id=ModelID)])))
  amb<-lookup[,.(n=uniqueN(model_id)),by=name][n>1L,name]
  lookup<-lookup[!name%in%amb]
  z<-fread(file.path(ROOT,"04_sDDS_model/SCP542_BC_counts.csv")); gs<-z[[1]]
  C<-as.matrix(z[,-1]);rm(z);gc()
  meta<-fread(file.path(ROOT,"04_sDDS_model/SCP542_BC_meta.csv"))
  stopifnot(setequal(colnames(C),meta$cell));meta<-meta[match(colnames(C),cell)]
  lines<-sort(unique(meta$cell_line))
  # Continuous linear expression; this is a cell-line mean, not integer UMI counts.
  bulk<-vapply(lines,function(s)rowMeans(C[,meta$cell_line==s,drop=FALSE]),numeric(nrow(C)))
  rownames(bulk)<-gs
  X1<-collapse_expression(t(bulk),map_symbols(gs,D$mapping),"linear")
  ncell<-as.integer(table(factor(meta$cell_line,levels=lines)))
  m1<-data.table(dataset="SCP542",sample_id=lines,cell_line=lines,n_cells=ncell,
       model_id=lookup$model_id[match(normalize(lines),lookup$name)])
  rm(C,bulk);gc()
  z<-fread(file.path(ROOT,"04_sDDS_model/GSE173634_pseudobulk_ensembl.csv"))
  gs<-z[[1]];X2<-t(as.matrix(z[,-1]));colnames(X2)<-gs
  X2<-collapse_expression(X2,map_symbols(gs,D$mapping),"linear")
  m2<-data.table(dataset="GSE173634",sample_id=rownames(X2),cell_line=rownames(X2),n_cells=NA_integer_,
       model_id=lookup$model_id[match(normalize(rownames(X2)),lookup$name)])
  mapping<-rbindlist(list(m1,m2)); mapping[,mapping_status:=ifelse(is.na(model_id),"unmapped_or_ambiguous","mapped")]
  write_table(mapping,"cross_assay_cell_line_mapping_v3.csv")
  obj<-list(SCP542=list(X=X1,meta=m1),GSE173634=list(X=X2,meta=m2))
  save_checkpoint(obj,"cross_assay_expression_v3.rds")
  obj
}
cross_assay<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  cp<-file.path(CHECK,"cross_assay_expression_v3.rds")
  assays<-if(file.exists(cp))readRDS(cp) else prepare_cross_assay(D)
  preds<-list()
  for(drug in D$eligibility$drug) {
    z<-readRDS(file.path(CHECK,paste0("trained_",drug,".rds")))
    idx<-which(D$labels$drug==drug);labs<-D$labels[idx]
    for(ds in names(assays)) {
      a<-assays[[ds]];R<-rank_samples(a$X,D$U)$ranks
      for(i in seq_len(nrow(a$meta))) {
        mid<-a$meta$model_id[i]; sid<-a$meta$sample_id[i]
        if(is.na(mid))next
        seen<-mid%in%z$final$bundles$xgboost_weighted$preprocessor$training_ids
        keys<-if(seen)names(z$outer)[vapply(z$outer,function(f)mid%in%f$test_ids,logical(1))] else "final"
        if(seen)stopifnot(length(keys)==5L)
        for(key in keys) {
          fits<-if(key=="final")z$final else z$outer[[key]]
          for(v in names(fits$bundles)) {
            b<-fits$bundles[[v]]
            stopifnot(!mid%in%b$preprocessor$training_ids)
            p<-predict_sdds(b,R[sid,,drop=FALSE],already_ranked=TRUE)
            p[,`:=`(dataset=ds,model_id=mid,cell_line=a$meta$cell_line[i],variant=v,outer_key=key,
              evaluation=if(seen)"cross_assay_heldout_cell_line" else "unseen_cell_line_final_model")]
            set(p,j="drug",value=drug)
            p[,y:=labs$y[match(model_id,labs$model_id)]]
            # Conflict/intermediate labels are not endpoints even if an assay was non-intermediate.
            status<-labs$label[match(mid,labs$model_id)]
            if(length(status)!=1L||is.na(status)||!status%in%c("sensitive","resistant"))p[,y:=NA_integer_]
            preds[[paste(drug,ds,sid,key,v)]]<-p
          }
          train_ids<-fits$bundles$xgboost_weighted$preprocessor$training_ids
          baseprob<-mean(labs$y[match(train_ids,labs$model_id)])
          base<-copy(p)
          base[,`:=`(variant="constant_baseline",probability_sensitive=baseprob)]
          preds[[paste(drug,ds,sid,key,"baseline")]]<-base
        }
      }
    }
  }
  out<-rbindlist(preds);write_table(out,"cross_assay_heldout_predictions_by_fold_v3.csv")
  agg<-out[,.(y=y[1],probability_sensitive=mean(probability_sensitive),
     coverage=min(coverage),n_heldout_models=.N,status=if(all(status=="accepted"))"accepted" else "rejected"),
      by=.(dataset,drug,model_id,cell_line,variant,evaluation)]
  write_table(agg,"cross_assay_cell_line_predictions_v3.csv")
  perf<-agg[!is.na(y)&status=="accepted",c(list(n=.N,n_sensitive=sum(y==1),n_resistant=sum(y==0)),
      as.list(performance(y,probability_sensitive))),by=.(dataset,drug,variant,evaluation)]
  write_table(perf,"cross_assay_cell_line_performance_v3.csv")
  print(perf[variant=="xgboost_weighted"])
  invisible(perf)
}
if(sys.nframe()==0L)cross_assay()
