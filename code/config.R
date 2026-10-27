options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
})
setDTthreads(2)
ROOT <- Sys.getenv("SPAGI_ROOT", file.path(getwd(),"data","SpaGI_BC"))
WORK <- Sys.getenv("SPAGI_WORK", file.path(getwd(),"work"))
OUT <- Sys.getenv("SPAGI_OUT", dirname(normalizePath(Sys.getenv("SPAGI_CODE_DIR","code"))))
CHECK <- Sys.getenv("SPAGI_CHECKPOINT", file.path(getwd(),"checkpoints"))
SEED <- 20261026L
VERSION <- "SpaGI-BC-v3-rank-20261026"
for (d in c(WORK,OUT,CHECK,file.path(OUT,c("tables","code","models","figures"))))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
sha <- function(f) digest::digest(file=f,algo="sha256",serialize=FALSE)
save_checkpoint <- function(x,name) {
  p <- file.path(WORK,name)
  saveRDS(x,p,compress="gzip")
  stopifnot(system2("cp",c(shQuote(p),shQuote(file.path(CHECK,name))))==0)
  invisible(p)
}
write_table <- function(x,name) fwrite(x,file.path(OUT,"tables",name),na="NA")
atlas_obs <- function(f=file.path(ROOT,"01_scRNA_processed/BC_scRNA_annotated.h5ad")) {
  decode <- function(n) {
    p<-paste0("obs/",n)
    cats<-as.character(h5read(f,paste0(p,"/categories"))); codes<-as.integer(h5read(f,paste0(p,"/codes")))
    cats[ifelse(codes>=0,codes+1L,NA_integer_)]
  }
  z<-data.table(cell_id=as.character(h5read(f,"obs/_index")))
  for(n in c("celltype","malignancy","patient","subtype")) z[,(n):=decode(n)]
  for(n in c("tme_myCAF","tme_SPP1_TAM","tme_CD8_Tex")) z[,(n):=h5read(f,paste0("obs/",n))]
  z
}
atlas_subset <- function(genes,cells=NULL,linear=FALSE,batch=500L,
                        f=file.path(ROOT,"01_scRNA_processed/BC_scRNA_annotated.h5ad")) {
  allg<-as.character(h5read(f,"raw/var/_index")); allc<-as.character(h5read(f,"obs/_index"))
  genes<-intersect(genes,allg); gi<-match(genes,allg)-1L
  if(is.null(cells)) cells<-seq_along(allc)
  stopifnot(all(cells %in% seq_along(allc)))
  ptr<-as.numeric(h5read(f,"raw/X/indptr"))
  res<-matrix(0,length(cells),length(genes),dimnames=list(allc[cells],genes))
  # contiguous raw CSR reads; target cells are selected only after reading the chunk.
  for(start in seq.int(1L,length(allc),by=batch)) {
    end<-min(start+batch-1L,length(allc)); selected<-which(cells>=start & cells<=end)
    if(!length(selected)) next
    lo<-ptr[start]+1; hi<-ptr[end+1L]
    if(hi<lo) next
    ids<-as.integer(h5read(f,"raw/X/indices",start=lo,count=hi-lo+1))
    vals<-as.numeric(h5read(f,"raw/X/data",start=lo,count=hi-lo+1))
    rows<-rep.int(seq.int(start,end),diff(ptr[start:(end+1L)]))
    cols<-match(ids,gi); use<-!is.na(cols) & rows %in% cells[selected]
    if(linear) vals<-expm1(vals)
    res[cbind(match(rows[use],cells),cols[use])]<-vals[use]
  }
  res
}
