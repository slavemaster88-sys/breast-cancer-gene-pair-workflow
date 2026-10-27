source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "train_models.R"))
auc_rank<-function(y,p) {
  n1<-sum(y==1);n0<-sum(y==0)
  if(!n1|!n0)return(NA_real_)
  (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
validate_models<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  checks<-list()
  for(drug in D$eligibility$drug) {
    name<-paste0("trained_",drug,".rds")
    stopifnot(file.exists(file.path(CHECK,name)),sha(file.path(CHECK,name))==sha(file.path(WORK,name)))
    z<-readRDS(file.path(CHECK,name))
    idx<-which(D$labels$drug==drug&D$labels$training_status=="eligible")
    labs<-D$labels[idx]
    stopifnot(all(z$oof$y==labs$y[match(z$oof$model_id,labs$model_id)]))
    for(key in names(z$outer)) {
      f<-z$outer[[key]]
      stopifnot(!length(intersect(f$train_ids,f$test_ids)),all(f$train_ids%in%labs$model_id),
                setequal(c(f$train_ids,f$test_ids),labs$model_id))
      for(b in f$bundles) {
        stopifnot(identical(b$drug,drug),setequal(b$preprocessor$training_ids,f$train_ids))
      }
      for(k in 1:3) {
        a<-strsplit(f$inner_audit$train_ids[k],";",fixed=TRUE)[[1]]
        v<-strsplit(f$inner_audit$validation_ids[k],";",fixed=TRUE)[[1]]
        stopifnot(!length(intersect(a,v)),setequal(c(a,v),f$train_ids))
      }
    }
    cc<-check_interface(z$final$bundles$xgboost_weighted,D$X[labs$model_id[1:4],,drop=FALSE],paste0("final_",drug))
    set(cc,j="drug",value=drug)
    checks[[drug]]<-cc
    finalpath<-file.path(OUT,"models",paste0("sDDS_",drug,"_v3.rds"))
    stopifnot(sha(finalpath)==sha(file.path(WORK,basename(finalpath))))
    invisible(readRDS(finalpath))
  }
  write_table(rbindlist(checks),"final_model_invariance_checks_v3.csv")
  oof<-fread(file.path(OUT,"tables/drug_OOF_all_v3.csv"))
  stopifnot(!anyDuplicated(oof[,.(drug,variant,repeat_id,model_id)]))
  p<-oof[,.(y=y[1],p=mean(probability_sensitive)),by=.(drug,variant,model_id)]
  comparisons<-list()
  for(drug in unique(p$drug)) {
    index<-which(p$drug==drug); ds<-p[index]
    a<-ds[variant=="xgboost_weighted",.(model_id,y,p_primary=p)]
    for(other in c("xgboost_unweighted","elastic_net_weighted","constant_baseline")) {
      b<-ds[ds$variant==other,.(model_id,p_other=p)]
      z<-merge(a,b,by="model_id")
      set.seed(SEED)
      delta<-auc_rank(z$y,z$p_primary)-auc_rank(z$y,z$p_other)
      n<-nrow(z)
      boot<-replicate(2000L,{i<-sample.int(n,n,replace=TRUE)
        auc_rank(z$y[i],z$p_primary[i])-auc_rank(z$y[i],z$p_other[i])})
      comparisons[[paste(drug,other)]]<-data.table(drug=drug,comparison=paste("xgboost_weighted minus",other),
         n=n,AUROC_meanOOF_difference=delta,
         conditional_CI_lo=quantile(boot,.025,na.rm=TRUE,names=FALSE),
         conditional_CI_hi=quantile(boot,.975,na.rm=TRUE,names=FALSE),
         max_probability_difference=max(abs(z$p_primary-z$p_other)))
    }
  }
  write_table(rbindlist(comparisons),"drug_paired_model_comparisons_v3.csv")
  cal<-p[,{
    lp<-qlogis(pmin(pmax(p,1e-6),1-1e-6))
    if(sd(lp)<1e-10) list(n=.N,intercept=NA_real_,slope=NA_real_,status="constant_prediction") else {
      f<-suppressWarnings(glm(y~lp,family=binomial()))
      list(n=.N,intercept=unname(coef(f)[1]),slope=unname(coef(f)[2]),status=if(f$converged)"descriptive_fit" else "nonconverged")
    }
  },by=.(drug,variant)]
  write_table(cal,"drug_OOF_calibration_v3.csv")
  cat("Verified 200 outer folds, inner isolation, final bundles, labels and prediction invariance.\n")
  invisible(rbindlist(comparisons))
}
if(sys.nframe()==0L)validate_models()
