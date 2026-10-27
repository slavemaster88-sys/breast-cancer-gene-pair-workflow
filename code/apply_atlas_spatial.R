source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "prepare_training.R"))
source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "train_models.R"))
load_primary_bundles<-function() {
  fs<-list.files(file.path(OUT,"models"),pattern="^sDDS_.*_v3.rds$",full.names=TRUE)
  z<-lapply(fs,function(f)readRDS(f)$bundles$xgboost_weighted)
  setNames(z,vapply(z,`[[`,character(1),"drug"))
}
apply_atlas<-function(max_batches=Inf) {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  obs<-atlas_obs(); ids<-which(obs$malignancy=="Malignant")
  bundles<-load_primary_bundles();sums<-list();signed<-list();counts<-list()
  for(drug in names(bundles)) {
    sums[[drug]]<-setNames(numeric(nrow(bundles[[drug]]$preprocessor$pairs)),bundles[[drug]]$preprocessor$feature_order)
    signed[[drug]]<-sums[[drug]];counts[[drug]]<-0L
  }
  predictions<-list(); bnum<-0L
  for(start in seq.int(1L,length(ids),by=250L)) {
    bnum<-bnum+1L;if(bnum>max_batches)break
    ii<-ids[start:min(start+249L,length(ids))]
    X<-atlas_subset(D$U,cells=ii,linear=FALSE)
    R<-rank_samples(X,D$U)$ranks
    for(drug in names(bundles)) {
      b<-bundles[[drug]]
      p<-predict_sdds(b,R,already_ranked=TRUE)
      p[,`:=`(patient=obs$patient[ii],subtype=obs$subtype[ii])]
      predictions[[paste(start,drug)]]<-p
      F<-transform_pair_features(R,b$preprocessor,already_ranked=TRUE)
      ok<-F$meta$status=="accepted"
      if(any(ok)) {
        phi<-as.matrix(predict(b$model,xgb.DMatrix(F$x[ok,,drop=FALSE]),predcontrib=TRUE))
        stopifnot(ncol(phi)==ncol(F$x)+1L)
        # Contributions are on model log-odds; bias excluded from feature rankings.
        phi<-phi[,seq_len(ncol(F$x)),drop=FALSE]
        sums[[drug]]<-sums[[drug]]+colSums(abs(phi))
        signed[[drug]]<-signed[[drug]]+colSums(phi);counts[[drug]]<-counts[[drug]]+nrow(phi)
      }
    }
  }
  pred<-rbindlist(predictions)
  imp<-rbindlist(lapply(names(bundles),function(drug)data.table(
    drug=drug,pair_id=names(sums[[drug]]),mean_abs_SHAP=sums[[drug]]/counts[[drug]],
    mean_SHAP=signed[[drug]]/counts[[drug]],n_cells=counts[[drug]],
    status="full_fit_descriptive_explanation")))
  suffix<-if(is.finite(max_batches))"_pilot" else ""
  write_table(pred,paste0("atlas_sDDS_cells_v3",suffix,".csv"))
  write_table(imp,paste0("atlas_sDDS_SHAP_summary_v3",suffix,".csv"))
  if(!is.finite(max_batches)) {
    means<-pred[,.(n_cells=.N,mean_sDDS=mean(probability_sensitive),coverage_min=min(coverage)),by=.(drug,patient,subtype)]
    write_table(means,"atlas_sDDS_patient_subtype_v3.csv")
    save_checkpoint(list(n=nrow(pred),SHAP=imp,patient_means=means),"atlas_sDDS_v3.rds")
  }
  cat("Atlas application:",length(unique(pred$sample_id)),"cells;",nrow(pred),"drug-cell scores\n")
  invisible(imp)
}
read_visium_counts<-function(h5,mapping=NULL) {
  genes<-as.character(h5read(h5,"matrix/features/name")); ids<-as.character(h5read(h5,"matrix/barcodes"))
  shp<-as.integer(h5read(h5,"matrix/shape"))
  C<-sparseMatrix(i=as.integer(h5read(h5,"matrix/indices")),
    p=as.integer(h5read(h5,"matrix/indptr")),x=as.numeric(h5read(h5,"matrix/data")),
    dims=shp,index1=FALSE,dimnames=list(genes,ids))
  if(!is.null(mapping)) {
    sy<-map_symbols(genes,mapping);ok<-which(!is.na(sy));ug<-sort(unique(sy[ok]))
    A<-sparseMatrix(i=match(sy[ok],ug),j=ok,x=1,dims=c(length(ug),length(genes)))
    C<-A%*%C;rownames(C)<-ug
  }
  C
}
moran_array6<-function(score,coords) {
  if(any(!is.finite(score))||sd(score)==0)return(NA_real_)
  kn<-RANN::nn2(coords,coords,k=7L)$nn.idx
  # Exclude own barcode rather than assuming the self neighbor has column 1.
  kn<-t(vapply(seq_len(nrow(kn)),function(i)head(kn[i,kn[i,]!=i],6L),integer(6)))
  z<-score-mean(score)
  sum(z*rowMeans(matrix(z[kn],nrow=nrow(kn))))/sum(z*z)
}
paired_v3<-function(low,high) {
  ok<-is.finite(low)&is.finite(high);d<-low[ok]-high[ok]
  exact<-!any(d==0)&!anyDuplicated(abs(d))
  if(!length(d))return(data.table(n=0L,median_difference=NA_real_,p=NA_real_,exact=FALSE))
  data.table(n=length(d),excluded=sum(!ok),n_lower=sum(d<0),n_higher=sum(d>0),n_zero=sum(d==0),
    median_difference=median(d),mean_difference=mean(d),
    p=if(all(d==0))1 else wilcox.test(d,mu=0,exact=exact,correct=FALSE)$p.value,exact=exact)
}
apply_spatial<-function(one_section=FALSE) {
  source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "spatial_pilot.R"))
  D<-readRDS(file.path(WORK,"training_data_v3.rds"));bundle<-load_primary_bundles()$Paclitaxel
  modules<-readRDS(file.path(WORK,"modules_v3.rds"))$modules
  module<-sort(unique(na.omit(modules$gene)))
  fs<-list.files(file.path(ROOT,"00_raw_data/Spatial_extra/extracted"),
                 pattern="filtered_feature_bc_matrix.h5$",full.names=TRUE)
  if(one_section)fs<-fs[1L]
  rows<-list();gates<-list()
  for(h5 in fs) {
    section<-strsplit(basename(h5),"_",fixed=TRUE)[[1]][2]
    legacy<-legacy_spatial_section(h5);frozen<-legacy$data
    old<-fread(file.path(ROOT,"06_clinical_validation",paste0("spot_ESV_",section,".csv")))
    setnames(old,1,"barcode");chk<-merge(frozen,old[,.(barcode,old_ESV=ESV)],by="barcode",all=TRUE)
    stopifnot(nrow(chk)==nrow(frozen),!anyNA(chk$old_ESV))
    pass<-all(abs(chk$legacy_ESV-chk$old_ESV)<=1e-5+1e-5*abs(chk$old_ESV))
    gates[[section]]<-data.table(section=section,n_spots=nrow(chk),max_abs_difference=max(abs(chk$legacy_ESV-chk$old_ESV)),pass=pass)
    if(!pass)stop("Spatial reproduction failed in ",section,"; publication paused")
    ad<-legacy$ad
    gs<-intersect(module,as.character(py_to_r(ad$var_names$to_numpy())))
    if(!length(gs))stop("No new module genes in ",section)
    sc$tl$score_genes(ad,r_to_py(as.list(gs)),score_name="ESV_module_v3",use_raw=FALSE)
    module_score<-py_to_r(ad$obs)[["ESV_module_v3"]]
    C<-read_visium_counts(h5,D$mapping)
    X<-t(as.matrix(C[match(D$U,rownames(C)),,drop=FALSE]));colnames(X)<-D$U
    p<-predict_sdds(bundle,X)
    stopifnot(identical(p$sample_id,frozen$barcode))
    sp<-copy(frozen)
    sp[,`:=`(section=section,module_score_v3=module_score,ESV_v3=module_score*TME_weight,
      sDDS_Paclitaxel_v3=p$probability_sensitive,coverage=p$coverage,prediction_status=p$status)]
    pospath<-sub("_filtered_feature_bc_matrix.h5$","_tissue_positions_list.csv.gz",h5)
    positions<-fread(pospath,header=FALSE)
    setnames(positions,c("barcode","in_tissue","array_row","array_col","pixel_row","pixel_col"))
    pos<-positions[match(sp$barcode,barcode)]
    stopifnot(!anyNA(pos$array_row),!anyDuplicated(pos$barcode))
    coords<-as.matrix(pos[,.(array_row,array_col)])
    sp[,`:=`(array_row=pos$array_row,array_col=pos$array_col)]
    write_table(sp,paste0("spatial_",section,"_scores_and_frozen_regions_v3.csv"))
    low<-sp$low_malignant;high<-sp$high_malignant
    rows[[section]]<-data.table(section=section,tissue=if(section%in%c("M1","M7","M9","M10"))"Normal" else "Cancer",
      n_spots=nrow(sp),n_low=sum(low),n_high=sum(high),
      ESV_low=mean(sp$ESV_v3[low]),ESV_high=mean(sp$ESV_v3[high]),
      sDDS_low=mean(sp$sDDS_Paclitaxel_v3[low]),sDDS_high=mean(sp$sDDS_Paclitaxel_v3[high]),
      Moran_ESV=moran_array6(sp$ESV_v3,coords),Moran_sDDS=moran_array6(sp$sDDS_Paclitaxel_v3,coords),
      module_n=length(gs),min_coverage=min(p$coverage))
    cat(section,"old-score reproduction passed;",nrow(sp),"spots scored.\n")
  }
  summary<-rbindlist(rows);gate<-rbindlist(gates)
  suffix<-if(one_section)"_pilot" else ""
  write_table(summary,paste0("spatial_section_summary_v3",suffix,".csv"))
  write_table(gate,paste0("spatial_all_reproduction_gates_v3",suffix,".csv"))
  if(!one_section) {
    tests<-rbindlist(lapply(c("ESV","sDDS"),function(v)cbind(endpoint=v,paired_v3(summary[[paste0(v,"_low")]],summary[[paste0(v,"_high")]]))))
    tests[,`:=`(q_BH=p.adjust(p,"BH"),unit="section; donor independence unconfirmed")]
    write_table(tests,"spatial_paired_statistics_v3.csv")
    save_checkpoint(list(summary=summary,gates=gate,statistics=tests),"spatial_v3.rds")
    print(tests)
  }
  invisible(summary)
}
if(sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if("--atlas"%in%args)apply_atlas(if("--pilot"%in%args)1L else Inf)
  if("--spatial"%in%args)apply_spatial("--pilot"%in%args)
}
