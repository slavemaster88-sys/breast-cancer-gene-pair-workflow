source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "clinical_reference.R"))
source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "train_models.R"))
source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "validate_models.R"))
suppressPackageStartupMessages(library(survival))
fit_clinical_transform<-function(d,features) {
  mu<-mean(d$module_mean_rank,na.rm=TRUE);sd0<-sd(d$module_mean_rank,na.rm=TRUE)
  stopifnot(is.finite(sd0),sd0>0)
  d<-copy(d);d[,ESV:=(module_mean_rank-mu)/sd0*TME_mean_proxy]
  X<-as.matrix(d[,..features]);storage.mode(X)<-"double"
  med<-apply(X,2,function(x)if(any(is.finite(x)))median(x[is.finite(x)]) else NA_real_)
  stopifnot(all(is.finite(med)))
  list(features=features,module_mean=mu,module_sd=sd0,median=med,training_ids=d$ResearchID,
       rule="fold-trained module mean/SD; training-fold predictor medians; glmnet internal scaling")
}
transform_clinical<-function(d,p) {
  d<-copy(d);d[,ESV:=(module_mean_rank-p$module_mean)/p$module_sd*TME_mean_proxy]
  X<-as.matrix(d[,p$features,with=FALSE]);storage.mode(X)<-"double"
  for(j in seq_len(ncol(X)))X[!is.finite(X[,j]),j]<-p$median[j]
  rownames(X)<-d$ResearchID;X
}
ridge_path<-function(X,y,lambda=NULL) {
  if(is.null(lambda)) {
    f<-glmnet(X,y,family="binomial",alpha=0,standardize=TRUE,nlambda=5)
    lambda<-max(f$lambda)*exp(seq(0,log(1e-4),length.out=50))
  }
  list(model=glmnet(X,y,family="binomial",alpha=0,standardize=TRUE,lambda=lambda),lambda=lambda)
}
fit_pcr<-function(d,features,seed) {
  folds<-stratified_folds(d$y,3L,seed); losses<-matrix(NA_real_,3,50)
  for(k in 1:3) {
    tr<-d[folds!=k];va<-d[folds==k]
    p<-fit_clinical_transform(tr,features)
    X<-transform_clinical(tr,p);V<-transform_clinical(va,p)
    fit<-ridge_path(X,tr$y)
    pr<-as.matrix(predict(fit$model,newx=V,s=fit$lambda,type="response"))
    losses[k,]<-apply(pr,2,function(z)loss_log(va$y,z))
  }
  mu<-colMeans(losses);se<-apply(losses,2,sd)/sqrt(3);best<-which.min(mu)
  idx<-which(mu<=mu[best]+se[best])[1]
  p<-fit_clinical_transform(d,features);X<-transform_clinical(d,p);fit<-ridge_path(X,d$y)
  list(preprocessor=p,model=fit$model,lambda=fit$lambda[idx],
    tuning=data.table(index=1:50,mean_logloss=mu,se=se,selected=seq_len(50)==idx))
}
predict_pcr<-function(fit,d)as.numeric(predict(fit$model,newx=transform_clinical(d,fit$preprocessor),s=fit$lambda,type="response"))
prepare_pcr<-function(C) {
  l<-parse_soft_v3(file.path(ROOT,"00_raw_data/Clinical/GSE194040_family.soft.gz"))
  l[,ResearchID:=sub("^ISPY2_","",title)]
  l<-l[,.(ResearchID,HR=as.numeric(hr),HER2=as.numeric(her2),MP=as.numeric(mp),y=as.numeric(pcr),arm)]
  a<-merge(l,C$features$ISPY194040,by.x="ResearchID",by.y="sample_id",all=FALSE)
  b<-fread(file.path(ROOT,"00_raw_data/Clinical/Pusztai2021_biomarkers.csv.gz"))
  b[,ResearchID:=as.character(ResearchID)]
  b<-b[,.(ResearchID,HR=as.numeric(get("HR.status..0.HR....1..HR..")),
      HER2=as.numeric(get("HER2.status..0..HER2...1..HER2..")),
      MP=as.numeric(get("MP.1.2..class..0..MP1..1..MP2.")),y=as.numeric(pCR.status),arm=Arm)]
  b<-merge(b,C$features$ISPY173839,by.x="ResearchID",by.y="sample_id",all=FALSE)
  a<-a[y%in%c(0,1)];b<-b[y%in%c(0,1)]
  stopifnot(!anyDuplicated(a$ResearchID),!anyDuplicated(b$ResearchID),!length(intersect(a$ResearchID,b$ResearchID)))
  stopifnot(all(a$MP%in%c(0,1)),all(b$MP%in%c(0,1)))
  write_table(data.table(variable="MP",included=TRUE,
       reason="Official GEO Overall design verified in both cohorts: 0=MP1 high risk; 1=MP2 (ultra) high risk"),
       "pCR_covariate_availability_v3.csv")
  write_table(data.table(cohort=c("GSE194040","GSE173839"),
    source_url=c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194040",
                 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173839"),
    coding="0=MP1; 1=MP2",verification_date="2026-10-26",status="compatible_source_documentation"),
    "MP_source_coding_verification_v3.csv")
  features<-list(M0=c("HR","HER2","MP"),M1=c("HR","HER2","MP","sDDS_Paclitaxel","sDDS_Olaparib"))
  features$M2<-c(features$M1,if(length(C$module))"ESV",if(length(C$tcell))"Tcell_SL")
  write_table(rbindlist(lapply(names(features),function(m)data.table(model=m,feature=features[[m]]))),"pCR_locked_model_features_v3.csv")
  list(train=a,test=b,features=features)
}
pcr_models<-function(pilot=FALSE) {
  C<-readRDS(file.path(WORK,"clinical_features_v3.rds"));P<-prepare_pcr(C)
  a<-P$train;b<-P$test;preds<-list();final<-list()
  reps<-if(pilot)1L else 5L;kmax<-if(pilot)1L else 5L
  for(rep in seq_len(reps)) {
    fold<-stratified_folds(a$y,5L,SEED+70000L+rep)
    for(k in seq_len(kmax)) {
      tr<-a[fold!=k];va<-a[fold==k]
      stopifnot(!length(intersect(tr$ResearchID,va$ResearchID)))
      for(m in names(P$features)) {
        fit<-fit_pcr(tr,P$features[[m]],SEED+70000L+rep*100L+k)
        p<-predict_pcr(fit,va)
        preds[[paste(rep,k,m)]]<-data.table(ResearchID=va$ResearchID,model=m,repeat_id=rep,fold=k,
          y=va$y,probability_pCR=p,training_n=nrow(tr),training_ids=paste(tr$ResearchID,collapse=";"))
      }
    }
  }
  out<-rbindlist(preds);suffix<-if(pilot)"_pilot" else ""
  fold_info<-unique(out[,.(repeat_id,fold,training_n,training_ids)])
  stopifnot(nrow(fold_info)==reps*kmax)
  write_table(fold_info,paste0("pCR_fold_training_IDs_v3",suffix,".csv"))
  write_table(out[,setdiff(names(out),"training_ids"),with=FALSE],paste0("pCR_OOF_predictions_v3",suffix,".csv"))
  perf<-out[,c(list(n=.N),as.list(performance(y,probability_pCR))),by=.(model,repeat_id)]
  write_table(perf,paste0("pCR_OOF_performance_by_repeat_v3",suffix,".csv"))
  if(pilot){print(perf);return(invisible(perf))}
  meanpred<-out[,.(y=y[1],p=mean(probability_pCR)),by=.(model,ResearchID)]
  summary<-perf[,lapply(.SD,mean),by=model,.SDcols=c("AUROC","AUPRC","Brier","logloss","prevalence")]
  cis<-meanpred[,{
    n<-.N;set.seed(SEED)
    boot<-replicate(2000L,{i<-sample.int(n,n,replace=TRUE);auc_rank(y[i],p[i])})
    .(n=n,AUROC_mean_OOF=auc_rank(y,p),
      conditional_CI_lo=quantile(boot,.025,names=FALSE,na.rm=TRUE),
      conditional_CI_hi=quantile(boot,.975,names=FALSE,na.rm=TRUE))
  },by=model]
  summary<-merge(summary,cis,by="model");write_table(summary,"pCR_OOF_summary_v3.csv")
  transfer<-list();coefs<-list()
  for(m in names(P$features)) {
    fit<-fit_pcr(a,P$features[[m]],SEED+79999L)
    final[[m]]<-fit
    pr<-predict_pcr(fit,b)
    transfer[[m]]<-data.table(ResearchID=b$ResearchID,model=m,y=b$y,arm=b$arm,probability_pCR=pr)
    cf<-as.matrix(coef(fit$model,s=fit$lambda))
    coefs[[m]]<-data.table(model=m,feature=rownames(cf),coefficient=cf[,1])
  }
  trans<-rbindlist(transfer)
  perftrans<-rbindlist(lapply(c("all","DOP"),function(s) {
    d<-if(s=="all")trans else trans[tolower(trimws(arm))=="durvalumab/olaparib"]
    d[,c(list(scope=s,n=.N,n_pCR=sum(y==1)),as.list(performance(y,probability_pCR))),by=model]
  }))
  write_table(trans,"pCR_frozen_transfer_predictions_v3.csv")
  write_table(perftrans,"pCR_frozen_transfer_performance_v3.csv")
  write_table(rbindlist(coefs),"pCR_final_model_coefficients_v3.csv")
  # Internal treatment-arm adjustment uses fixed indicator levels from the training source.
  armlevels<-sort(unique(a$arm));B<-model.matrix(~factor(arm,levels=armlevels)-1,data=a)
  colnames(B)<-make.names(colnames(B))
  armdata<-cbind(a,as.data.table(B));rownames(armdata)<-NULL
  sensitivity<-list()
  for(m in names(P$features)) {
    for(rep in 1:5) {
      fold<-stratified_folds(a$y,5L,SEED+70000L+rep)
      for(k in 1:5) {
        fit<-fit_pcr(armdata[fold!=k],c(P$features[[m]],colnames(B)),SEED+70000L+rep*100L+k)
        va<-armdata[fold==k]
        sensitivity[[paste(m,rep,k)]]<-data.table(ResearchID=va$ResearchID,model=m,repeat_id=rep,
          fold=k,y=va$y,probability_pCR=predict_pcr(fit,va))
      }
    }
  }
  ss<-rbindlist(sensitivity);write_table(ss,"pCR_arm_adjusted_internal_OOF_v3.csv")
  write_table(ss[,as.list(performance(y,probability_pCR)),by=.(model,repeat_id)],
              "pCR_arm_adjusted_internal_performance_v3.csv")
  save_checkpoint(list(final=final,features=P$features,training_ids=a$ResearchID,
     module_training_mean=final$M2$preprocessor$module_mean,module_training_sd=final$M2$preprocessor$module_sd,
     OOF=out,transfer=trans,summary=summary),"clinical_pCR_models_v3.rds")
  print(summary);print(perftrans)
  invisible(summary)
}
survival_models<-function() {
  C<-readRDS(file.path(WORK,"clinical_features_v3.rds"))
  PCR<-readRDS(file.path(WORK,"clinical_pCR_models_v3.rds"))
  results<-list();allph<-list();allcoefs<-list();allcases<-list()
  for(tag in c("METABRIC","TCGA")) {
    fn<-if(tag=="METABRIC")"METABRIC_clinical.csv" else "TCGA_PANCAN_clinical.csv"
    cl<-fread(file.path(ROOT,"00_raw_data/Clinical",fn));setnames(cl,1,"patient_id")
    manifest<-copy(C$sample_manifests[[tag]])
    manifest[,selection_status:=ifelse(sampleType=="Primary Solid Tumor","candidate","non_primary")]
    setorder(manifest,patientId,sampleId)
    manifest[selection_status=="candidate",selection_status:=ifelse(seq_len(.N)==1L,"selected_primary","duplicate_primary_excluded"),by=patientId]
    write_table(manifest,paste0("survival_",tag,"_sample_selection_v3.csv"))
    selected<-manifest[selection_status=="selected_primary",.(patient_id=patientId,sample_id=sampleId)]
    stopifnot(!anyDuplicated(selected$patient_id))
    # Existing patient-aggregated expression is valid only where the source has one sample.
    multi<-manifest[,.(n=.N),by=patientId][n>1L,patientId]
    if(any(selected$patient_id%in%multi))stop("Multiple source samples per patient require sample-resolved expression retrieval")
    features<-copy(C$features[[tag]]);setnames(features,"sample_id","patient_id")
    d<-merge(selected,cl,by="patient_id",all.x=TRUE)
    d<-merge(d,features,by="patient_id",all.x=TRUE)
    d[,ESV_raw:=(module_mean_rank-PCR$module_training_mean)/PCR$module_training_sd*TME_mean_proxy]
    if(tag=="METABRIC") {
      d[,`:=`(age=as.numeric(AGE_AT_DIAGNOSIS),stage=as.character(TUMOR_STAGE),subtype=as.character(CLAUDIN_SUBTYPE))]
    } else {
      d[,`:=`(age=as.numeric(AGE),stage=sub("^STAGE ([IV]+)[ABC]?$","\\1",toupper(trimws(AJCC_PATHOLOGIC_TUMOR_STAGE))),
        subtype=as.character(SUBTYPE))]
      d[,stage:=as.character(match(stage,c("I","II","III","IV")))]
    }
    d[!stage%in%as.character(1:4),stage:=NA_character_]
    d[is.na(subtype)|subtype%in%c("","NA","Not Available","[Not Available]"),subtype:=NA_character_]
    endpoints<-if(tag=="TCGA")c("OS","DSS") else "OS"
    for(ep in endpoints) {
      z<-copy(d)
      z[,`:=`(time=suppressWarnings(as.numeric(get(paste0(ep,"_MONTHS")))),
        event=ifelse(grepl("^1",get(paste0(ep,"_STATUS"))),1L,
                     ifelse(grepl("^0",get(paste0(ep,"_STATUS"))),0L,NA_integer_)))]
      vars<-c("time","event","ESV_raw","age","stage","subtype","myCAF_proxy","SPP1_TAM_proxy")
      z[,exclusion_reason:=fcase(!is.finite(time)|time<=0,"missing_or_nonpositive_time",
          is.na(event),"missing_event",!complete.cases(z[,vars,with=FALSE]),"incomplete_prespecified_covariates",default="included")]
      write_table(z[,.(patient_id,sample_id,time,event,age,stage,subtype,ESV_raw,myCAF_proxy,SPP1_TAM_proxy,exclusion_reason)],
                  paste0("survival_",tag,"_",ep,"_case_flow_v3.csv"))
      z<-z[exclusion_reason=="included"];z[,`:=`(stage=factor(stage,levels=as.character(1:4)),subtype=factor(subtype))]
      z[,ESV_z:=(ESV_raw-mean(ESV_raw))/sd(ESV_raw)]
      sdESV<-sd(z$ESV_raw)
      full_formula<-Surv(time,event)~ESV_z+age+stage+subtype+myCAF_proxy+SPP1_TAM_proxy
      fit<-coxph(full_formula,data=z,x=TRUE,y=TRUE)
      ph<-cox.zph(fit)
      pht<-data.table(term=rownames(ph$table),as.data.table(ph$table))
      pht[,`:=`(cohort=tag,endpoint=ep)]
      allph[[paste(tag,ep)]]<-pht
      violate<-intersect(c("stage","subtype"),pht[p<.05,term])
      esv_time<-pht[term=="ESV_z",p]<.05
      terms<-c("ESV_z","age",if("stage"%in%violate)"strata(stage)" else "stage",
        if("subtype"%in%violate)"strata(subtype)" else "subtype","myCAF_proxy","SPP1_TAM_proxy")
      if(esv_time)terms<-c(terms,"tt(ESV_z)")
      form<-as.formula(paste("Surv(time,event)~",paste(terms,collapse="+")))
      ttfunc<-function(x,t,...)x*log(pmax(t,1e-6))
      adjusted<-coxph(form,data=z,tt=ttfunc,x=TRUE,y=TRUE)
      reduced<-coxph(as.formula(paste("Surv(time,event)~",paste(terms[!grepl("ESV_z",terms,fixed=TRUE)],collapse="+"))),
                    data=z)
      cf<-coef(adjusted);V<-vcov(adjusted);dfESV<-if(esv_time)2L else 1L
      lrt<-max(0,2*(adjusted$loglik[2]-reduced$loglik[2]));p<-pchisq(lrt,dfESV,lower.tail=FALSE)
      epv<-sum(z$event)/length(cf)
      times<-if(esv_time)c(12,60) else NA_real_
      for(tm in times) {
        vec<-setNames(numeric(length(cf)),names(cf));vec["ESV_z"]<-1
        if(esv_time)vec["tt(ESV_z)"]<-log(tm)
        est<-sum(vec*cf);se<-sqrt(as.numeric(t(vec)%*%V%*%vec))
        results[[paste(tag,ep,tm)]]<-data.table(cohort=tag,endpoint=ep,n=nrow(z),events=sum(z$event),
          coefficient_df=length(cf),events_per_df=epv,adequate_EPV=epv>=10,
          PH_ESV_time_dependent=esv_time,stratified_terms=paste(violate,collapse=";"),time_months=tm,
          HR_per_SD=exp(est),CI_lo=exp(est-1.96*se),CI_hi=exp(est+1.96*se),
          HR_per_raw_unit=exp(est/sdESV),raw_CI_lo=exp((est-1.96*se)/sdESV),raw_CI_hi=exp((est+1.96*se)/sdESV),
          ESV_cohort_SD=sdESV,overall_ESV_p=p,
          interpretation=if(epv>=10)"exploratory_adjusted_association" else "low_EPV_no_confirmatory_adjusted_P")
      }
      tab<-summary(adjusted)$coefficients
      allcoefs[[paste(tag,ep)]]<-data.table(cohort=tag,endpoint=ep,term=rownames(tab),as.data.table(tab))
      allcases[[paste(tag,ep)]]<-z[,.(patient_id,sample_id,time,event,ESV_raw,ESV_z,age,stage,subtype,myCAF_proxy,SPP1_TAM_proxy)]
    }
  }
  res<-rbindlist(results)
  tests<-unique(res[,.(cohort,endpoint,overall_ESV_p)])
  tests[,overall_ESV_q_BH:=p.adjust(overall_ESV_p,"BH")]
  res<-merge(res,tests,by=c("cohort","endpoint","overall_ESV_p"),all.x=TRUE)
  write_table(res,"survival_ESV_results_v3.csv");write_table(rbindlist(allph),"survival_Schoenfeld_diagnostics_v3.csv")
  write_table(rbindlist(allcoefs),"survival_adjusted_coefficients_v3.csv")
  write_table(rbindlist(allcases,idcol="cohort_endpoint"),"survival_complete_case_data_v3.csv")
  save_checkpoint(list(results=res,PH=rbindlist(allph),coefs=rbindlist(allcoefs)),"survival_v3.rds")
  print(res);invisible(res)
}
if(sys.nframe()==0L) {
  source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "validate_models.R"))
  if("--pcr"%in%commandArgs(TRUE))pcr_models("--pilot"%in%commandArgs(TRUE))
  if("--survival"%in%commandArgs(TRUE))survival_models()
}
