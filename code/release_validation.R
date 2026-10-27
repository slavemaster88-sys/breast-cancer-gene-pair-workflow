# Packaging-only tests using real existing samples, never simulated data.
Sys.setenv(SPAGI_CODE_DIR=Sys.getenv("SPAGI_CODE_DIR",file.path(getwd(),"code")))
source(file.path(Sys.getenv("SPAGI_CODE_DIR"),"sdds_interface.R"))
validate_release<-function(){
  scripts<-list.files(Sys.getenv("SPAGI_CODE_DIR"),pattern="\\.R$",full.names=TRUE)
  for(f in scripts)invisible(parse(f))
  D<-readRDS(file.path(WORK,"training_data_v3.rds"));X<-D$X[1:3,,drop=FALSE]
  checks<-list()
  for(drug in unique(D$eligibility$drug)){
    f<-file.path(OUT,"models",paste0("sDDS_",drug,"_v3.rds"))
    B<-readRDS(f)$bundles$xgboost_weighted
    p<-predict_sdds(B,X)
    q<-predict_sdds(B,X[1,,drop=FALSE])
    orig<-readRDS(file.path(WORK,basename(f)))$bundles$xgboost_weighted
    base<-predict_sdds(orig,X)
    delta<-max(abs(p$probability_sensitive-base$probability_sensitive))
    batch<-abs(q$probability_sensitive-p$probability_sensitive[1])
    stopifnot(delta==0,batch==0,all(p$status=="accepted"))
    checks[[drug]]<-data.table(test="portable_shared_interface_real_samples",scope=drug,n=nrow(X),
        max_prediction_difference=delta,max_single_batch_difference=batch,pass=TRUE)
  }
  source(file.path(Sys.getenv("SPAGI_CODE_DIR"),"clinical_models.R"))
  C<-readRDS(file.path(WORK,"clinical_features_v3.rds"));P<-prepare_pcr(C)
  F<-readRDS(file.path(OUT,"models/pCR_frozen_models_v3.rds"))
  obs<-fread(file.path(OUT,"tables/pCR_frozen_transfer_predictions_v3.csv"))
  for(m in names(F$final)){
    p<-predict_pcr(F$final[[m]],P$test)
    original<-obs[model==m][match(P$test$ResearchID,ResearchID)]
    delta<-max(abs(p-original$probability_pCR))
    stopifnot(delta<1e-14)
    checks[[m]]<-data.table(test="frozen_clinical_bundle_transfer_readback",scope=m,n=length(p),
      max_prediction_difference=delta,max_single_batch_difference=NA_real_,pass=TRUE)
  }
  # Independent AUROC implementation confirms primary means from saved predictions.
  oof<-fread(file.path(OUT,"tables/drug_OOF_all_v3.csv"))[variant=="xgboost_weighted"]
  recalc<-oof[,.(auc=as.numeric(pROC::auc(pROC::roc(y,probability_sensitive,levels=c(0,1),direction="<",quiet=TRUE)))),by=.(drug,repeat_id)]
  perf<-fread(file.path(OUT,"tables/drug_performance_by_repeat_v3.csv"))[variant=="xgboost_weighted"]
  z<-merge(recalc,perf,by=c("drug","repeat_id"))
  stopifnot(nrow(z)==40L,max(abs(z$auc-z$AUROC))<1e-14)
  checks[["AUROC"]]<-data.table(test="independent_pROC_vs_saved_rank_AUC",scope="eight_drug_40_repeats",
     n=40,max_prediction_difference=max(abs(z$auc-z$AUROC)),max_single_batch_difference=NA_real_,pass=TRUE)
  # Verify DOCX XML and embedding count without reading staging copies.
  docs<-list.files(OUT,pattern="\\.docx$",full.names=TRUE)
  for(f in docs){
    td<-file.path(WORK,paste0("inspect_",sub("\\.docx$","",basename(f))))
    if(dir.exists(td))unlink(td,recursive=TRUE);dir.create(td)
    utils::unzip(f,exdir=td)
    x<-xml2::read_xml(file.path(td,"word/document.xml"));ns<-xml2::xml_ns(x)
    text<-paste(xml2::xml_text(xml2::xml_find_all(x,"//w:t",ns)),collapse=" ")
    ndraw<-length(xml2::xml_find_all(x,"//w:drawing",ns))
    stopifnot(nchar(text)>500)
    if(grepl("Embedded",f))stopifnot(ndraw==6L)
    checks[[basename(f)]]<-data.table(test="DOCX_XML_and_content_readback",scope=basename(f),
       n=ndraw,max_prediction_difference=NA_real_,max_single_batch_difference=NA_real_,pass=TRUE)
  }
  tests<-rbindlist(checks,use.names=TRUE)
  fwrite(tests,file.path(OUT,"tables/release_packaging_checks_v3.csv"))
  print(tests)
}
if(sys.nframe()==0L)validate_release()
