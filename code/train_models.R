source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "sdds_interface.R"))
suppressPackageStartupMessages({library(glmnet);library(xgboost)})
loss_log<-function(y,p) {
  p<-pmin(pmax(p,1e-8),1-1e-8)
  -mean(y*log(p)+(1-y)*log1p(-p))
}
stratified_folds<-function(y,k,seed) {
  set.seed(seed); out<-integer(length(y))
  for(v in sort(unique(y))) {
    ix<-which(y==v); out[sample(ix,length(ix))]<-rep(seq_len(k),length.out=length(ix))
  }
  out
}
xgb_fit<-function(X,y,depth,nrounds,seed) {
  xgb.train(params=list(objective="binary:logistic",eval_metric="logloss",
     eta=.05,subsample=.8,colsample_bytree=.8,min_child_weight=3,lambda=1,
     max_depth=as.integer(depth),nthread=2L,seed=as.integer(seed)),
     data=xgb.DMatrix(X,label=y),nrounds=as.integer(nrounds),verbose=0)
}
enet_path<-function(X,y,lambda=NULL) {
  if(is.null(lambda)) {
    base<-glmnet(X,y,family="binomial",alpha=.5,standardize=TRUE,nlambda=5)
    lambda<-max(base$lambda)*exp(seq(0,log(1e-4),length.out=50L))
  }
  list(model=glmnet(X,y,family="binomial",alpha=.5,standardize=TRUE,lambda=lambda,
                   thresh=1e-9,maxit=100000L),lambda=lambda)
}
make_bundle<-function(pre,model,drug,family,weighted,selected,lambda=NULL) {
  list(version=VERSION,drug=drug,family=family,weighted=weighted,preprocessor=pre,
       model=model,lambda=lambda,selected_parameters=selected,seed=SEED,
       software=c(R=as.character(getRversion()),xgboost=as.character(packageVersion("xgboost")),
                  glmnet=as.character(packageVersion("glmnet"))))
}
fit_outer_models<-function(D,ids,drug,seed) {
  # Avoid data.table's parent-name ambiguities in programmatic selection.
  index<-which(D$labels$drug==drug & D$labels$training_status=="eligible")
  labs<-D$labels[index]
  stopifnot(all(labs$drug==drug),!anyDuplicated(labs$model_id))
  y<-labs$y[match(ids,labs$model_id)]; stopifnot(!anyNA(y),!anyDuplicated(ids))
  R<-D$ranks[ids,,drop=FALSE]
  fold<-stratified_folds(y,3L,seed)
  grid<-CJ(depth=1:3,nrounds=c(100L,200L))
  xloss<-matrix(NA_real_,3L,nrow(grid)); eloss<-matrix(NA_real_,3L,50L)
  inner_audit<-list()
  for(k in 1:3) {
    tr<-fold!=k;va<-fold==k
    pre<-fit_preprocessor(R[tr,,drop=FALSE],D$U,D$pairs,already_ranked=TRUE)
    train<-transform_pair_features(R[tr,,drop=FALSE],pre,already_ranked=TRUE)$x
    test<-transform_pair_features(R[va,,drop=FALSE],pre,already_ranked=TRUE)$x
    stopifnot(length(intersect(rownames(train),rownames(test)))==0)
    for(j in seq_len(nrow(grid))) {
      fit<-xgb_fit(train,y[tr],grid$depth[j],grid$nrounds[j],seed+k)
      p<-as.numeric(predict(fit,xgb.DMatrix(test)))
      xloss[k,j]<-loss_log(y[va],p)
    }
    ep<-enet_path(train,y[tr])
    pp<-as.matrix(predict(ep$model,newx=test,s=ep$lambda,type="response"))
    eloss[k,]<-apply(pp,2,function(p)loss_log(y[va],p))
    inner_audit[[k]]<-data.table(inner_fold=k,train_ids=paste(ids[tr],collapse=";"),
      validation_ids=paste(ids[va],collapse=";"),n_features=ncol(train),preprocessor_hash=pre$spec_hash)
  }
  grid[,mean_logloss:=colMeans(xloss)]
  best<-grid[order(mean_logloss,depth,nrounds)][1]
  em<-colMeans(eloss); es<-apply(eloss,2,sd)/sqrt(3)
  emin<-which.min(em); eidx<-which(em<=em[emin]+es[emin])[1L]
  pre<-fit_preprocessor(R,D$U,D$pairs,already_ranked=TRUE)
  bundles<-list()
  primary_lambda<-NULL
  for(weighted in c(TRUE,FALSE)) {
    tag<-if(weighted)"weighted" else "unweighted"
    X<-transform_pair_features(R,pre,weighted=weighted,already_ranked=TRUE)$x
    xf<-xgb_fit(X,y,best$depth,best$nrounds,seed)
    bundles[[paste0("xgboost_",tag)]]<-make_bundle(pre,xf,drug,"xgboost",weighted,as.list(best))
    ep<-enet_path(X,y,lambda=primary_lambda)
    if(weighted) primary_lambda<-ep$lambda
    bundles[[paste0("elastic_net_",tag)]]<-make_bundle(pre,ep$model,drug,"elastic_net",weighted,
       list(alpha=.5,lambda_rule="inner-fold-local-50-point-lambda.1se",
            lambda_index=eidx,lambda_fraction=exp(seq(0,log(1e-4),length.out=50))[eidx]),
       lambda=primary_lambda[eidx])
  }
  list(bundles=bundles,xgb_grid=grid,enet_tuning=data.table(index=1:50,mean_logloss=em,se=es,selected=seq_len(50)==eidx),
       inner_audit=rbindlist(inner_audit))
}
check_interface<-function(bundle,X,tag="pilot") {
  stopifnot(nrow(X)>=3L)
  p<-predict_sdds(bundle,X)
  one<-predict_sdds(bundle,X[1,,drop=FALSE])
  reordered<-predict_sdds(bundle,X[rev(seq_len(nrow(X))),rev(seq_len(ncol(X))),drop=FALSE])
  reload_path<-file.path(WORK,paste0(tag,"_reload_test.rds"));saveRDS(bundle,reload_path)
  reload<-predict_sdds(readRDS(reload_path),X)
  altered<-X; altered[-1,]<-0;alteredp<-predict_sdds(bundle,altered)
  diffs<-c(single_batch=abs(one$probability_sensitive-p$probability_sensitive[1]),
    row_column_order=max(abs(reordered$probability_sensitive[match(p$sample_id,reordered$sample_id)]-p$probability_sensitive)),
    reload=max(abs(reload$probability_sensitive-p$probability_sensitive)),
    target_batch_independence=abs(alteredp$probability_sensitive[1]-p$probability_sensitive[1]))
  # Rejected low coverage is explicit, never a zero score.
  low<-predict_sdds(bundle,X[,seq_len(floor(.9*ncol(X))),drop=FALSE])
  stopifnot(all(is.na(low$probability_sensitive)),all(low$status=="coverage_below_95pct"))
  stopifnot(all(is.finite(diffs)),all(diffs<1e-12))
  data.table(test=names(diffs),max_abs_difference=as.numeric(diffs),pass=diffs<1e-12)
}
pilot_model<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  d<-D$labels[drug=="Paclitaxel" & training_status=="eligible"]
  f<-stratified_folds(d$y,5L,SEED)
  tr<-d$model_id[f!=1L];te<-d$model_id[f==1L]
  t<-proc.time()
  fitted<-fit_outer_models(D,tr,"Paclitaxel",SEED+100L)
  elapsed<-unname((proc.time()-t)[3])
  predictions<-rbindlist(lapply(names(fitted$bundles),function(n) {
    b<-fitted$bundles[[n]]; p<-predict_sdds(b,D$X[te,,drop=FALSE]); p[,variant:=n]; p
  }))
  predictions[,y:=d$y[match(sample_id,d$model_id)]]
  checks<-rbindlist(lapply(names(fitted$bundles),function(n)check_interface(fitted$bundles[[n]],D$X[te,,drop=FALSE],n)[,variant:=n]))
  # XGB uses order-only splits; standardized EN can exactly cancel fixed positive weights.
  write_table(predictions,"pilot_Paclitaxel_predictions.csv")
  write_table(checks,"pilot_interface_invariance_checks.csv")
  timings<-data.table(drug="Paclitaxel",outer_fold=1L,n_train=length(tr),n_test=length(te),
    n_features=nrow(fitted$bundles$xgboost_weighted$preprocessor$pairs),elapsed_seconds=elapsed,
    projected_208_fits_hours=elapsed*(8*25+8)/3600)
  mem<-grep("^VmHWM:",readLines("/proc/self/status"),value=TRUE)
  timings[,process_peak_RSS_MB:=as.numeric(gsub("[^0-9]","",mem))/1024]
  timings[,R_heap_highwater_MB:=sum(gc()[,6])]
  timings[,memory_measurement_note:="R heap high-water excludes native-library memory; process RSS unavailable when NA"]
  write_table(timings,"pilot_model_runtime.csv")
  save_checkpoint(list(fitted=fitted,train_ids=tr,test_ids=te,predictions=predictions,timings=timings),"pilot_model_v3.rds")
  print(timings);print(checks)
  invisible(timings)
}
performance<-function(y,p) {
  ok<-is.finite(y)&is.finite(p); y<-y[ok];p<-p[ok]
  if(length(unique(y))<2)return(c(AUROC=NA,AUPRC=NA,
     Brier=if(length(y))mean((p-y)^2) else NA_real_,
     logloss=if(length(y))loss_log(y,p) else NA_real_,prevalence=mean(y)))
  auc<-as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))
  d<-data.table(y=y,p=p)[,.(pos=sum(y),n=.N),by=p][order(-p)]
  ap<-sum((d$pos/sum(y))*(cumsum(d$pos)/cumsum(d$n)))
  c(AUROC=auc,AUPRC=ap,Brier=mean((p-y)^2),logloss=loss_log(y,p),prevalence=mean(y))
}
train_all<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  allpred<-list(); runtimes<-list()
  for(di in seq_along(D$eligibility$drug)) {
    drug<-D$eligibility$drug[di]
    if(D$eligibility$status[di]!="evaluable")next
    done<-file.path(CHECK,paste0("trained_",drug,".rds"))
    if(file.exists(done)) {
      z<-readRDS(done);allpred[[drug]]<-z$oof;cat("Reused verified completed stage:",drug,"\n");next
    }
    clock<-proc.time()[3]
    index<-which(D$labels$drug==drug & D$labels$training_status=="eligible")
    labs<-D$labels[index]
    stopifnot(all(labs$drug==drug),!anyDuplicated(labs$model_id))
    preds<-list();outer<-list(); audits<-list()
    kfold<-D$eligibility$outer_folds[di]
    for(rep in 1:5) {
      fold<-stratified_folds(labs$y,kfold,SEED+di*1000L+rep)
      for(k in seq_len(kfold)) {
        key<-paste0("repeat",rep,"_fold",k)
        trainids<-labs$model_id[fold!=k]; testids<-labs$model_id[fold==k]
        cache<-file.path(CHECK,paste0("fold_",drug,"_",key,".rds"))
        if(file.exists(cache)) {
          fit<-readRDS(cache)
          stopifnot(identical(fit$train_ids,trainids),identical(fit$test_ids,testids))
        } else {
          fit<-fit_outer_models(D,trainids,drug,SEED+di*10000L+rep*100L+k)
          fit$train_ids<-trainids;fit$test_ids<-testids
          save_checkpoint(fit,basename(cache))
          stopifnot(sha(cache)==sha(file.path(WORK,basename(cache))))
        }
        stopifnot(!length(intersect(trainids,testids)))
        for(n in names(fit$bundles)) {
          p<-predict_sdds(fit$bundles[[n]],D$ranks[testids,,drop=FALSE],already_ranked=TRUE)
          p[,`:=`(variant=n,model_id=sample_id,repeat_id=rep,fold=k,
            y=labs$y[match(sample_id,labs$model_id)],train_n=length(trainids),
            training_ids=paste(trainids,collapse=";"))]
          preds[[paste(key,n)]]<-p
        }
        base<-data.table(sample_id=testids,coverage=1,n_missing_genes=0,status="accepted",drug=drug,
          version=VERSION,probability_sensitive=mean(labs$y[fold!=k]),variant="constant_baseline",
          model_id=testids,repeat_id=rep,fold=k,y=labs$y[fold==k],train_n=length(trainids),
          training_ids=paste(trainids,collapse=";"))
        preds[[paste(key,"baseline")]]<-base
        outer[[key]]<-fit
        cat(drug,key,"completed at",format(Sys.time()),"\n")
      }
    }
    final<-fit_outer_models(D,labs$model_id,drug,SEED+di*10000L+999L)
    index<-which(D$cutoffs$drug==drug)
    final$label_definition<-D$cutoffs[index]
    for(n in names(final$bundles)) {
      final$bundles[[n]]$gene_mapping<-D$mapping
      final$bundles[[n]]$input_manifest<-D$input_manifest
      final$bundles[[n]]$label_definition<-final$label_definition
    }
    saveRDS(final,file.path(WORK,paste0("sDDS_",drug,"_v3.rds")),compress="gzip")
    stopifnot(system2("cp",c(shQuote(file.path(WORK,paste0("sDDS_",drug,"_v3.rds"))),
                           shQuote(file.path(OUT,"models",paste0("sDDS_",drug,"_v3.rds")))))==0)
    z<-list(drug=drug,oof=rbindlist(preds),outer=outer,final=final,elapsed_seconds=proc.time()[3]-clock)
    save_checkpoint(z,paste0("trained_",drug,".rds"))
    stopifnot(sha(done)==sha(file.path(WORK,basename(done))))
    write_table(z$oof,paste0("OOF_",drug,"_v3.csv"))
    allpred[[drug]]<-z$oof
    runtimes[[drug]]<-data.table(drug=drug,elapsed_seconds=z$elapsed_seconds)
  }
  oof<-rbindlist(allpred);write_table(oof,"drug_OOF_all_v3.csv")
  perf<-oof[,as.list(performance(y,probability_sensitive)),by=.(drug,variant,repeat_id)]
  write_table(perf,"drug_performance_by_repeat_v3.csv")
  means<-perf[,lapply(.SD,mean),by=.(drug,variant),.SDcols=c("AUROC","AUPRC","Brier","logloss","prevalence")]
  meanpred<-oof[,.(y=y[1],probability_sensitive=mean(probability_sensitive)),by=.(drug,variant,model_id)]
  ci<-meanpred[,{
    set.seed(SEED)
    boot<-replicate(2000L,{
      ix<-sample(seq_len(.N),.N,replace=TRUE)
      performance(y[ix],probability_sensitive[ix])["AUROC"]
    })
    .(n=.N,AUROC_mean_OOF=performance(y,probability_sensitive)["AUROC"],
      AUROC_conditional_lo=quantile(boot,.025,na.rm=TRUE,names=FALSE),
      AUROC_conditional_hi=quantile(boot,.975,na.rm=TRUE,names=FALSE))
  },by=.(drug,variant)]
  summary<-merge(means,ci,by=c("drug","variant"))
  write_table(summary,"drug_performance_summary_v3.csv")
  write_table(rbindlist(runtimes),"drug_training_runtime_v3.csv")
  print(summary)
  invisible(summary)
}
if(sys.nframe()==0L) {
  if("--all"%in%commandArgs(TRUE))train_all() else pilot_model()
}
