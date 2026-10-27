source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
Sys.setenv(RETICULATE_PYTHON=Sys.getenv("SPAGI_PYTHON",Sys.which("python")))
suppressPackageStartupMessages(library(reticulate))
sc<-import("scanpy",convert=FALSE)
legacy_spatial_section<-function(h5) {
  ad<-sc$read_10x_h5(h5)
  ad$var_names_make_unique()
  sc$pp$normalize_total(ad,target_sum=1e4)
  sc$pp$log1p(ad)
  genes<-as.character(py_to_r(ad$var_names$to_numpy()))
  sigs<-jsonlite::fromJSON(file.path(ROOT,"05_spatial_niche/signature_genes.json"))
  oldtri<-fread(file.path(ROOT,"03_TME_extrinsic_SV/Tripartite_LR_SV_network.csv"))
  oldmod<-sort(unique(c(oldtri$receptor,unlist(strsplit(na.omit(oldtri$receptor_upstream_SV_hubs),";",fixed=TRUE)))))
  for(n in c("Malignant","myCAF","SPP1_TAM","CD8_T")) {
    gs<-intersect(sigs[[n]],genes);stopifnot(length(gs)>0)
    sc$tl$score_genes(ad,r_to_py(as.list(gs)),score_name=paste0("s_",n),use_raw=FALSE)
  }
  sc$tl$score_genes(ad,r_to_py(as.list(intersect(oldmod,genes))),score_name="esv_mod",use_raw=FALSE)
  o<-py_to_r(ad$obs); ids<-rownames(o)
  tw<-(o$s_myCAF+o$s_SPP1_TAM)/2
  tw<-(tw-min(tw))/(max(tw)-min(tw)+1e-9)
  q<-quantile(o$s_Malignant,c(.25,.75),names=FALSE)
  d<-data.table(barcode=ids,malignant_score=o$s_Malignant,legacy_ESV=o$esv_mod*tw,
      low_malignant=o$s_Malignant<=q[1],high_malignant=o$s_Malignant>=q[2],TME_weight=tw)
  list(data=d,ad=ad)
}
spatial_pilot<-function() {
  f<-list.files(file.path(ROOT,"00_raw_data/Spatial_extra/extracted"),pattern="_M1_filtered_feature_bc_matrix.h5$",full.names=TRUE)
  stopifnot(length(f)==1L)
  clock<-proc.time()[3]
  x<-legacy_spatial_section(f);d<-x$data
  old<-fread(file.path(ROOT,"06_clinical_validation/spot_ESV_M1.csv"))
  setnames(old,1,"barcode")
  chk<-merge(d,old[,.(barcode,old_ESV=ESV)],by="barcode",all=TRUE)
  stopifnot(nrow(chk)==nrow(d),!anyNA(chk$old_ESV))
  absdiff<-max(abs(chk$legacy_ESV-chk$old_ESV))
  passed<-all(abs(chk$legacy_ESV-chk$old_ESV)<=1e-5+1e-5*abs(chk$old_ESV))
  summary<-fread(file.path(ROOT,"06_clinical_validation/spatial_replication_gse213688.csv"))[sample=="M1"]
  low<-mean(d[low_malignant==TRUE,legacy_ESV]);high<-mean(d[high_malignant==TRUE,legacy_ESV])
  meanpass<-abs(low-summary$ESV_boundary)<=1e-5+1e-5*abs(summary$ESV_boundary) &
             abs(high-summary$ESV_core)<=1e-5+1e-5*abs(summary$ESV_core)
  gate<-data.table(section="M1",n_spots=nrow(d),max_spot_abs_difference=absdiff,
    legacy_low_mean=low,saved_low_mean=summary$ESV_boundary,
    legacy_high_mean=high,saved_high_mean=summary$ESV_core,
    spot_pass=passed,mean_pass=meanpass,elapsed_seconds=proc.time()[3]-clock)
  write_table(gate,"pilot_spatial_reproduction_gate.csv")
  print(gate)
  if(!passed|!meanpass) stop("Legacy spatial reproduction gate failed; publication paused")
  write_table(d,"pilot_M1_frozen_regions.csv")
  save_checkpoint(list(regions=d,gate=gate),"pilot_spatial_v3.rds")
  invisible(gate)
}
if(sys.nframe()==0L)spatial_pilot()
