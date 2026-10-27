source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "sdds_interface.R"))
DRUGS<-c("Paclitaxel","Olaparib","Talazoparib","Palbociclib","Doxorubicin","Afatinib","Lapatinib","Bortezomib")
build_mapping<-function() {
  f<-gzfile(Sys.getenv("SPAGI_GENCODE",file.path(getwd(),"data","reference","gencode.v44.annotation.gtf.gz")),"rt")
  on.exit(close(f)); parts<-list(); i<-0L
  repeat {
    lines<-readLines(f,n=100000L)
    if(!length(lines)) break
    lines<-lines[grepl("\tgene\t",lines,fixed=TRUE)]
    if(!length(lines)) next
    i<-i+1L
    parts[[i]]<-data.table(
      ensembl=sub('.*gene_id "([^"]+)".*','\\1',lines),
      symbol=sub('.*gene_name "([^"]+)".*','\\1',lines))
  }
  m<-unique(rbindlist(parts)); m[,ensembl:=sub("\\..*$","",ensembl)]
  m[,ensembl_ambiguous:=uniqueN(symbol)>1L,by=ensembl]
  m[,symbol_multiple_ensembl:=uniqueN(ensembl)>1L,by=symbol]
  write_table(m,"GENCODE_v44_symbol_mapping.csv")
  m
}
map_symbols<-function(x,m) {
  x<-trimws(as.character(x)); clean<-sub("\\..*$","",x)
  ens<-grepl("^ENSG[0-9]+",clean)
  out<-x; out[ens]<-m$symbol[match(clean[ens],m$ensembl)]
  amb<-m[ensembl_ambiguous==TRUE,ensembl]
  out[ens & clean%in%amb]<-NA_character_
  out[!out%in%m$symbol]<-NA_character_
  out
}
collapse_expression<-function(X,symbols,scale=c("log","linear")) {
  scale<-match.arg(scale); X<-as.matrix(X)
  keep<-which(!is.na(symbols)&nzchar(symbols)); X<-X[,keep,drop=FALSE]; sy<-symbols[keep]
  if(!anyDuplicated(sy)) {colnames(X)<-sy;return(X)}
  groups<-split(seq_along(sy),sy)
  out<-vapply(groups,function(j) {
    if(length(j)==1L) X[,j] else if(scale=="linear") rowSums(X[,j,drop=FALSE],na.rm=FALSE)
    else apply(X[,j,drop=FALSE],1,median,na.rm=TRUE)
  },numeric(nrow(X)))
  rownames(out)<-rownames(X); out
}
prepare_training<-function() {
  m<-build_mapping()
  cclepath<-file.path(ROOT,"00_raw_data/DepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
  h<-names(fread(cclepath,nrows=0))
  avail<-list(CCLE=sub(" \\([^)]+\\)$","",h[-1]),
    atlas=h5read(file.path(ROOT,"01_scRNA_processed/BC_scRNA_annotated.h5ad"),"raw/var/_index"),
    I_SPY194040=fread(file.path(ROOT,"00_raw_data/Clinical/ISPY2_geneLevel.txt.gz"),skip=1,header=FALSE,select=1)[[1]],
    I_SPY173839=fread(file.path(ROOT,"00_raw_data/Clinical/Pusztai2021_geneLevel.txt.gz"),skip=1,header=FALSE,select=1)[[1]],
    SCP542=fread(file.path(ROOT,"04_sDDS_model/SCP542_BC_counts.csv"),select=1)[[1]],
    GSE173634=fread(file.path(ROOT,"04_sDDS_model/GSE173634_pseudobulk_ensembl.csv"),select=1)[[1]])
  fs<-list.files(file.path(ROOT,"00_raw_data/Spatial_extra/extracted"),pattern="filtered_feature_bc_matrix.h5$",full.names=TRUE)
  for(f in fs) avail[[sub("_filtered.*","",basename(f))]]<-h5read(f,"matrix/features/name")
  mapping_audit<-rbindlist(lapply(names(avail),function(ds) {
    x<-avail[[ds]]; mapped<-map_symbols(x,m)
    data.table(dataset=ds,input_identifier=x,symbol=mapped,status=ifelse(is.na(mapped),"unmapped_or_ambiguous",
       ifelse(duplicated(mapped)|duplicated(mapped,fromLast=TRUE),"duplicate_symbol","mapped")))
  }))
  write_table(mapping_audit,"assay_gene_mapping_audit.csv")
  available<-lapply(avail,function(x) unique(na.omit(map_symbols(x,m))))
  U<-sort(Reduce(intersect,available))
  stopifnot(length(U)>=1000L)
  write_table(data.table(gene=U),"sDDS_common_gene_universe.csv")
  write_table(data.table(assay=names(avail),n_identifiers=lengths(avail),n_mapped_symbols=lengths(available),
                         common_U=length(U)),"assay_gene_availability.csv")
  pairs<-readRDS(file.path(WORK,"network_v3.rds"))$pairs
  pairs[,universe_status:=fifelse(symbol1%in%U & symbol2%in%U & is.finite(W)&W>0,"eligible","outside_U_or_invalid_W")]
  write_table(pairs,"sDDS_fixed_pair_panel.csv")
  pairs<-pairs[universe_status=="eligible"];stopifnot(nrow(pairs)>0)
  models<-fread(file.path(ROOT,"00_raw_data/DepMap/Model.csv"))
  breast<-models[OncotreeLineage=="Breast"]
  D<-fread(cclepath)
  D<-D[D[[1]]%in%breast$ModelID]; ids<-D[[1]]
  X<-as.matrix(D[,-1]); rownames(X)<-ids
  X<-collapse_expression(X,map_symbols(avail$CCLE,m),"log")[,U,drop=FALSE]
  ranks<-rank_samples(X,U)$ranks
  # Re-extract reference records rather than reusing old label decisions.
  lookup<-rbindlist(list(models[,.(alias=toupper(CellLineName),model_id=ModelID)],
                         models[,.(alias=toupper(StrippedCellLineName),model_id=ModelID)],
                         models[,.(alias=toupper(ModelID),model_id=ModelID)]))
  lookup<-unique(lookup[!is.na(alias)&nzchar(alias)])
  amb<-lookup[,.(n=uniqueN(model_id)),by=alias][n>1,alias]
  lookup<-lookup[!alias%in%amb]
  response<-rbindlist(lapply(c("GDSC1","GDSC2"),function(ds) {
    d<-fread(file.path(ROOT,"00_raw_data/DrugResponse",paste0(ds,"_fitted.csv")))
    d<-d[tolower(DRUG_NAME)%in%tolower(DRUGS)]
    d[,model_id:=lookup$model_id[match(toupper(CELL_LINE_NAME),lookup$alias)]]
    d<-d[model_id%in%breast$ModelID]
    d[,.(dataset=ds,drug=DRUGS[match(tolower(DRUG_NAME),tolower(DRUGS))],
      model_id,cell_line=CELL_LINE_NAME,value=AUC,metric="AUC",source_record=seq_len(.N))]
  }))
  ti<-fread(file.path(ROOT,"00_raw_data/DrugResponse/PRISM_treatment_info.csv"))
  lfcpath<-file.path(ROOT,"00_raw_data/DrugResponse/PRISM_primary_LFC.csv")
  lh<-names(fread(lfcpath,nrows=0))
  ti<-ti[tolower(name)%in%tolower(DRUGS) & column_name%in%lh]
  lfc<-fread(lfcpath,select=c(lh[1],unique(ti$column_name)))
  cli<-fread(file.path(ROOT,"00_raw_data/DrugResponse/PRISM_cell_line_info.csv"))
  prismids<-lfc[[1]]
  alt<-!startsWith(prismids,"ACH"); prismids[alt]<-cli$depmap_id[match(prismids[alt],cli$row_name)]
  response<-rbindlist(c(list(response),lapply(DRUGS,function(drug) {
    cols<-ti[tolower(name)==tolower(drug),column_name]
    if(!length(cols)) return(NULL)
    vals<-rowMeans(as.matrix(lfc[,..cols]),na.rm=TRUE)
    data.table(dataset="PRISM19Q4",drug=drug,model_id=prismids,cell_line=lfc[[1]],value=vals,
      metric="mean_LFC_across_existing_doses_screens",source_record=seq_along(vals))[
      model_id%in%breast$ModelID & is.finite(value)]
  })),use.names=TRUE)
  response<-response[is.finite(value)]
  # Retain every experimental response record in the reference definition and trace.
  cuts<-response[,.(q1=quantile(value,1/3,names=FALSE,type=7),
      q2=quantile(value,2/3,names=FALSE,type=7),n_reference_records=.N,n_reference_models=uniqueN(model_id)),by=.(dataset,drug)]
  response<-merge(response,cuts,by=c("dataset","drug"),all.x=TRUE,sort=FALSE)
  response[,label:=fcase(value<q1,"sensitive",value>q2,"resistant",default="intermediate")]
  labels<-response[,{
    labs<-unique(label[label!="intermediate"])
    .(label=if(length(labs)==1L)labs else if(length(labs)>1L)"conflict" else "intermediate",
      n_source_records=.N,datasets=paste(sort(unique(dataset)),collapse=";"),
      source_labels=paste(paste(dataset,label,sep=":"),collapse=";"))
  },by=.(drug,model_id)]
  labels[,expression_available:=model_id%in%rownames(X)]
  labels[,training_status:=fcase(label=="conflict","cross_assay_conflict",label=="intermediate","intermediate",
    !expression_available,"no_training_expression",default="eligible")]
  labels[,y:=fifelse(label=="sensitive",1L,fifelse(label=="resistant",0L,NA_integer_))]
  labels[,cell_line:=models$CellLineName[match(model_id,models$ModelID)]]
  write_table(response,"breast_response_source_labels_v3.csv")
  write_table(cuts,"breast_response_cutoffs_v3.csv")
  write_table(labels,"breast_consensus_labels_v3.csv")
  counts<-labels[training_status=="eligible",.(n=.N,n_sensitive=sum(y==1L),n_resistant=sum(y==0L)),by=drug]
  counts[,outer_folds:=ifelse(pmin(n_sensitive,n_resistant)>=10L,5L,
                             ifelse(pmin(n_sensitive,n_resistant)>=6L,3L,NA_integer_))]
  counts[,status:=ifelse(is.na(outer_folds),"not_evaluable","evaluable")]
  write_table(counts,"drug_training_eligibility_v3.csv")
  input_manifest<-fread(file.path(OUT,"input_manifest_v3.csv"))
  obj<-list(U=U,pairs=pairs,mapping=m,X=X,ranks=ranks,labels=labels,eligibility=counts,
            cutoffs=cuts,input_manifest=input_manifest,seed=SEED,version=VERSION)
  save_checkpoint(obj,"training_data_v3.rds")
  cat("Common U:",length(U),"pairs:",nrow(pairs),"breast CCLE:",nrow(X),"\n");print(counts)
  invisible(obj)
}
if(sys.nframe()==0L) prepare_training()
