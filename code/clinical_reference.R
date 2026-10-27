source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "prepare_training.R"))
suppressPackageStartupMessages(library(httr2))
parse_soft_v3<-function(path) {
  x<-readLines(gzfile(path));out<-list();row<-list()
  for(l in x) {
    if(startsWith(l,"^SAMPLE = ")) {
      if(length(row))out[[length(out)+1L]]<-row
      row<-list(GSM=sub("^\\^SAMPLE = ","",l))
    }
    if(startsWith(l,"!Sample_title = "))row$title<-sub("^!Sample_title = ","",l)
    if(startsWith(l,"!Sample_characteristics_ch1 = ")) {
      v<-sub("^!Sample_characteristics_ch1 = ","",l)
      if(grepl(":",v,fixed=TRUE))row[[sub(":.*$","",v)]]<-sub("^[^:]+: *","",v)
    }
  }
  if(length(row))out[[length(out)+1L]]<-row
  rbindlist(out,fill=TRUE)
}
read_microarray<-function(name,mapping) {
  f<-file.path(ROOT,"00_raw_data/Clinical",name)
  hdr<-strsplit(readLines(gzfile(f),n=1),"\t",fixed=TRUE)[[1]]
  cols<-if(name=="ISPY2_geneLevel.txt.gz")c("gene",hdr) else c("gene",hdr[-1])
  d<-fread(f,skip=1,header=FALSE,col.names=cols)
  g<-d[[1]];X<-t(as.matrix(d[,-1]));colnames(X)<-g
  rownames(X)<-sub("^ISPY2_","",trimws(rownames(X)))
  stopifnot(!anyDuplicated(rownames(X)))
  collapse_expression(X,map_symbols(g,mapping),"log")
}
read_survival_expression<-function(tag,mapping) {
  names<-if(tag=="METABRIC")c("METABRIC_expression_hv.csv","METABRIC_expression_sdds_genes.csv","METABRIC_expression_tme.csv")
         else "TCGA_PANCAN_expression.csv"
  mats<-lapply(names,function(fn) {
    d<-fread(file.path(ROOT,"00_raw_data/Clinical",fn));g<-d[[1]]
    x<-t(as.matrix(d[,-1]));colnames(x)<-map_symbols(g,mapping);x
  })
  common<-Reduce(intersect,lapply(mats,rownames))
  X<-do.call(cbind,lapply(mats,function(x)x[common,,drop=FALSE]))
  # Repeated downloads of identical symbols are merged without target-cohort fitting.
  collapse_expression(X,colnames(X),"log")
}
cbio_get<-function(endpoint) {
  request(paste0("https://www.cbioportal.org/api/",endpoint)) |>
    req_timeout(90) |> req_retry(max_tries=2,max_seconds=150) |> req_perform() |>
    resp_body_json(simplifyVector=TRUE)
}
get_sample_manifest<-function(study) {
  p<-file.path(WORK,paste0(study,"_sample_manifest.rds"))
  if(file.exists(p))d<-readRDS(p) else {
    d<-cbio_get(paste0("studies/",study,"/samples?projection=DETAILED&pageSize=100000"))
    saveRDS(d,p)
  }
  d<-as.data.table(d);stopifnot(all(c("sampleId","patientId","sampleType")%in%names(d)))
  write_table(d,paste0(study,"_sample_manifest_v3.csv"))
  save_checkpoint(d,paste0(study,"_sample_manifest_v3.rds"));d
}
fetch_missing_genes<-function(tag,X,genes) {
  missing<-setdiff(genes,colnames(X))
  if(!length(missing))return(X)
  study<-if(tag=="METABRIC")"brca_metabric" else "brca_tcga_pan_can_atlas_2018"
  profile<-if(tag=="METABRIC")"brca_metabric_mrna" else "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
  cache<-file.path(CHECK,paste0(tag,"_new_module_gene_records_v3.rds"))
  if(file.exists(cache)) records<-readRDS(cache) else {
    hdr<-names(fread(file.path(ROOT,"00_raw_data/DepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv"),nrows=0))[-1]
    known<-data.table(symbol=sub(" \\([^)]+\\)$","",hdr),
                      entrezGeneId=as.integer(sub(".*\\(([0-9]+)\\)$","\\1",hdr)))
    ids<-lapply(missing,function(g)tryCatch({
      hit<-known[which(known$symbol==g)]
      if(nrow(hit)==1L && is.finite(hit$entrezGeneId))return(hit)
      r<-cbio_get(paste0("genes/",utils::URLencode(g,reserved=TRUE)))
      if(is.list(r) && !is.null(r$entrezGeneId))data.table(symbol=g,entrezGeneId=as.integer(r$entrezGeneId)) else NULL
    },error=function(e)NULL))
    ids<-rbindlist(ids)
    records<-data.table()
    if(nrow(ids)) {
      response<-request(paste0("https://www.cbioportal.org/api/molecular-profiles/",profile,"/molecular-data/fetch?projection=DETAILED")) |>
        req_body_json(list(sampleListId=paste0(study,"_all"),entrezGeneIds=as.list(ids$entrezGeneId))) |>
        req_timeout(150) |> req_retry(max_tries=2,max_seconds=240) |> req_perform() |>
        resp_body_json(simplifyVector=TRUE)
      records<-as.data.table(response)[,intersect(c("sampleId","patientId","entrezGeneId","value","molecularProfileId"),names(response)),with=FALSE]
      if(nrow(records))records[,symbol:=ids$symbol[match(entrezGeneId,ids$entrezGeneId)]]
    }
    save_checkpoint(records,basename(cache))
  }
  write_table(records,paste0(tag,"_retrieved_new_module_expression_v3.csv"))
  if(!nrow(records))return(X)
  records[,patient_key:=if(tag=="METABRIC")sampleId else substr(sampleId,1,12)]
  if(anyDuplicated(records[,.(symbol,patient_key)]))stop("New expression records have multiple samples per patient; sample-resolved selection required")
  wide<-dcast(records,patient_key~symbol,value.var="value")
  M<-as.matrix(wide[,-1]);rownames(M)<-wide$patient_key
  for(g in intersect(missing,colnames(M))) {
    X<-cbind(X,M[match(rownames(X),rownames(M)),g])
    colnames(X)[ncol(X)]<-g
  }
  X
}
reference_labels<-function() {
  f<-file.path(ROOT,"01_scRNA_processed/BC_scRNA_annotated.h5ad")
  ct<-as.character(h5read(f,"obs/w_celltype_minor/categories"))[
    as.integer(h5read(f,"obs/w_celltype_minor/codes"))+1L]
  obs<-atlas_obs();spp1<-obs$tme_SPP1_TAM
  mac<-ct%in%c("Macrophage","Monocyte","Cycling_Myeloid")
  thr<-quantile(spp1[mac],.75,na.rm=TRUE,names=FALSE)
  ct[mac&spp1>=thr]<-"SPP1_TAM";ct[mac&spp1<thr]<-"Macrophage_other"
  collapse<-c("Cancer LumA SC"="Malignant","Cancer LumB SC"="Malignant","Cancer Basal SC"="Malignant",
    "Cancer Her2 SC"="Malignant","Cancer Cycling"="Malignant","Myoepithelial"="Malignant",
    "Luminal Progenitors"="Normal_Epi","Mature Luminal"="Normal_Epi",
    "T cells CD4+"="CD4_T","T cells CD8+"="CD8_T","Cycling T-cells"="CD8_T","NKT cells"="CD8_T",
    "NK cells"="NK","B cells Memory"="B_cells","B cells Naive"="B_cells","Plasmablasts"="Plasma",
    "CAFs myCAF-like"="myCAF","CAFs MSC iCAF-like"="iCAF",
    "Endothelial ACKR1"="Endo","Endothelial CXCL12"="Endo","Endothelial RGS5"="Endo",
    "Endothelial Lymphatic LYVE1"="Endo","PVL Differentiated"="PVL","PVL Immature"="PVL","DCs"="DC")
  idx<-match(ct,names(collapse));ct[!is.na(idx)]<-unname(collapse[idx[!is.na(idx)]])
  obs[,reference_type:=ct]
  obs
}
rank_ref<-function(X,U) {
  stopifnot(all(U%in%colnames(X)))
  M<-X[,U,drop=FALSE]
  R<-t(apply(M,1,function(x) {
    n<-sum(is.finite(x));x[!is.finite(x)]<-NA
    (rank(x,ties.method="average",na.last="keep")-.5)/n
  }))
  dimnames(R)<-dimnames(M);R
}
composition_proxy<-function(R,basis) {
  stopifnot(identical(colnames(R),rownames(basis)))
  coefs<-matrix(NA_real_,nrow(R),ncol(basis),dimnames=list(rownames(R),colnames(basis)))
  for(i in seq_len(nrow(R))) {
    ok<-is.finite(R[i,])&apply(basis,1,function(x)all(is.finite(x)))
    if(sum(ok)<.95*nrow(basis))next
    f<-nnls::nnls(basis[ok,,drop=FALSE],R[i,ok])
    v<-as.numeric(coef(f));if(sum(v)>0)coefs[i,]<-v/sum(v)
  }
  coefs
}
prepare_clinical<-function() {
  D<-readRDS(file.path(WORK,"training_data_v3.rds"))
  X<-list(ISPY194040=read_microarray("ISPY2_geneLevel.txt.gz",D$mapping),
    ISPY173839=read_microarray("Pusztai2021_geneLevel.txt.gz",D$mapping),
    METABRIC=read_survival_expression("METABRIC",D$mapping),
    TCGA=read_survival_expression("TCGA",D$mapping))
  mods<-readRDS(file.path(WORK,"modules_v3.rds"))$modules
  module<-sort(unique(na.omit(mods$gene)))
  tcell<-readRDS(file.path(WORK,"network_v3.rds"))$tcell[selected==TRUE,gene]
  extra<-unique(c(module,tcell))
  missing<-rbindlist(lapply(names(X),function(tag)data.table(cohort=tag,gene=setdiff(extra,colnames(X[[tag]])))))
  write_table(missing,"clinical_new_gene_missing_before_retrieval.csv")
  for(tag in c("METABRIC","TCGA"))X[[tag]]<-fetch_missing_genes(tag,X[[tag]],extra)
  old<-fread(file.path(ROOT,"06_clinical_validation/deconv_reference_profiles.csv"))[[1]]
  common<-Reduce(intersect,lapply(X,colnames))
  Uref<-sort(unique(c(intersect(old,common),extra[extra%in%common])))
  all_extra<-extra%in%common
  availability<-data.table(gene=extra,in_module=extra%in%module,in_Tcell=extra%in%tcell,
     common_measurable=all_extra)
  write_table(availability,"clinical_new_gene_availability_v3.csv")
  # Missing requested genes are not silently replaced with partial modules.
  if(any(!all_extra))stop("New module/T-cell genes unavailable after targeted retrieval; clinical scoring paused")
  stopifnot(length(Uref)>=100L)
  write_table(data.table(gene=Uref),"clinical_reference_gene_universe_v3.csv")
  obs<-reference_labels();E<-atlas_subset(Uref,linear=TRUE)
  stopifnot(identical(rownames(E),as.character(obs$cell_id)))
  group<-interaction(obs$patient,obs$reference_type,drop=TRUE,lex.order=TRUE)
  sums<-rowsum(E,group,reorder=TRUE)
  nc<-table(group);means<-sums/as.numeric(nc[rownames(sums)])
  # Equal weighting of atlas patient labels within each frozen reference type.
  group_info<-unique(data.table(group_key=as.character(group),patient=obs$patient,type=obs$reference_type))
  types<-sort(unique(obs$reference_type));ref<-vapply(types,function(t) {
    keys<-group_info[type==t,group_key]
    colMeans(means[match(keys,rownames(means)),,drop=FALSE])
  },numeric(length(Uref)))
  rownames(ref)<-Uref
  basis<-t(rank_ref(t(ref),Uref))
  write_table(data.table(gene=rownames(basis),as.data.table(basis)),"clinical_rank_reference_basis_v3.csv")
  write_table(group_info,"reference_patient_type_contributors_v3.csv")
  features<-list()
  for(tag in names(X)) {
    R<-rank_ref(X[[tag]],Uref);P<-composition_proxy(R,basis)
    out<-data.table(sample_id=rownames(R),module_mean_rank=rowMeans(R[,module,drop=FALSE]),
      Tcell_SL=if(length(tcell))rowMeans(R[,tcell,drop=FALSE]) else NA_real_,
      myCAF_proxy=P[,"myCAF"],SPP1_TAM_proxy=P[,"SPP1_TAM"])
    out[,TME_mean_proxy:=(myCAF_proxy+SPP1_TAM_proxy)/2]
    write_table(out,paste0("clinical_raw_features_",tag,"_v3.csv"))
    write_table(data.table(sample_id=rownames(P),as.data.table(P)),paste0("clinical_composition_proxies_",tag,"_v3.csv"))
    features[[tag]]<-out
  }
  # Drug predictions use exactly the shared CCLE-trained interface, not clinical-cohort z scores.
  for(tag in c("ISPY194040","ISPY173839")) {
    for(drug in c("Paclitaxel","Olaparib")) {
      b<-readRDS(file.path(OUT,"models",paste0("sDDS_",drug,"_v3.rds")))$bundles$xgboost_weighted
      p<-predict_sdds(b,X[[tag]])
      stopifnot(identical(p$sample_id,features[[tag]]$sample_id),all(p$status=="accepted"))
      set(features[[tag]],j=paste0("sDDS_",drug),value=p$probability_sensitive)
    }
    write_table(features[[tag]],paste0("clinical_raw_features_",tag,"_v3.csv"))
  }
  manifest<-list(TCGA=get_sample_manifest("brca_tcga_pan_can_atlas_2018"),
                 METABRIC=get_sample_manifest("brca_metabric"))
  obj<-list(Uref=Uref,basis=basis,module=module,tcell=tcell,features=features,sample_manifests=manifest,
     reference_rule="raw linear atlas patient-label balanced mean per type, then per-profile Uref ranks",
     software=as.character(packageVersion("nnls")))
  save_checkpoint(obj,"clinical_features_v3.rds")
  cat("Clinical Uref:",length(Uref),"module:",length(module),"T-cell:",length(tcell),"\n")
  invisible(obj)
}
if(sys.nframe()==0L)prepare_clinical()
