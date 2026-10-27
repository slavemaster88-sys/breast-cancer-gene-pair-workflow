source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
rank_samples <- function(expression,U,min_coverage=.95) {
  X<-as.matrix(expression)
  stopifnot(!is.null(colnames(X)),!anyDuplicated(colnames(X)),!is.null(rownames(X)),!anyDuplicated(rownames(X)))
  pos<-match(U,colnames(X)); R<-matrix(NA_real_,nrow(X),length(U),dimnames=list(rownames(X),U))
  found<-which(!is.na(pos)); R[,found]<-X[,pos[found],drop=FALSE]
  coverage<-rowMeans(is.finite(R))
  R[!is.finite(R)]<-NA_real_
  for(i in seq_len(nrow(R))) {
    n<-sum(!is.na(R[i,]))
    if(n>0) R[i,]<-(rank(R[i,],ties.method="average",na.last="keep")-.5)/n
  }
  list(ranks=R,meta=data.table(sample_id=rownames(X),coverage=coverage,
       n_missing_genes=rowSums(is.na(R)),status=ifelse(coverage>=min_coverage,"accepted","coverage_below_95pct")))
}
fit_preprocessor <- function(expression,U,pairs,already_ranked=FALSE) {
  stopifnot(length(U)>=1000,!anyDuplicated(U),all(c("symbol1","symbol2","W","pair_id")%in%names(pairs)))
  if(already_ranked) {
    R<-as.matrix(expression)[,U,drop=FALSE]
    meta<-data.table(coverage=rowMeans(is.finite(R)))
  } else {
    q<-rank_samples(expression,U); R<-q$ranks; meta<-q$meta
  }
  stopifnot(all(meta$coverage>=.95))
  meds<-apply(R,2,function(x)if(any(is.finite(x))) median(x[is.finite(x)]) else NA_real_)
  for(j in seq_len(ncol(R))) R[!is.finite(R[,j]),j]<-meds[j]
  mu<-colMeans(R); sd0<-apply(R,2,sd)
  valid<-is.finite(sd0)&sd0>sqrt(.Machine$double.eps)
  pairs<-copy(as.data.table(pairs))
  pairs[,fold_exclusion:=fcase(!symbol1%in%U|!symbol2%in%U,"outside_U",
    !is.finite(W)|W<=0,"nonpositive_or_nonfinite_W",default="")]
  idx<-which(pairs$fold_exclusion=="" &
        !(pairs$symbol1%in%names(valid)[valid]&pairs$symbol2%in%names(valid)[valid]))
  pairs[idx,fold_exclusion:="zero_variance_or_unmeasured_gene"]
  active<-pairs[fold_exclusion==""]
  stopifnot(nrow(active)>0,!anyDuplicated(active$pair_id))
  pre<-list(version=VERSION,U=U,mean=mu,sd=sd0,median_rank=meds,pairs=active,
    pair_audit=pairs,min_coverage=.95,rank_rule="(average_tied_rank-.5)/measured_U_count",
    missing_rule="training_fold_gene_median_rank",training_ids=rownames(R))
  z<-sweep(sweep(R,2,mu,"-"),2,sd0,"/")
  F<-(z[,active$symbol1,drop=FALSE]+z[,active$symbol2,drop=FALSE])/2
  colnames(F)<-active$pair_id
  pv<-apply(F,2,sd); good<-is.finite(pv)&pv>sqrt(.Machine$double.eps)
  if(any(!good)) {
    pre$pair_audit[pair_id%in%active$pair_id[!good],fold_exclusion:="zero_variance_pair"]
    pre$pairs<-active[good]
  }
  stopifnot(nrow(pre$pairs)>0)
  pre$feature_order<-pre$pairs$pair_id
  pre$spec_hash<-digest::digest(pre,algo="sha256")
  pre
}
transform_pair_features <- function(expression,preprocessor,weighted=TRUE,already_ranked=FALSE) {
  p<-preprocessor
  if(already_ranked) {
    R<-as.matrix(expression)[,p$U,drop=FALSE]
    meta<-data.table(sample_id=rownames(R),coverage=rowMeans(is.finite(R)),n_missing_genes=rowSums(!is.finite(R)))
    meta[,status:=ifelse(coverage>=p$min_coverage,"accepted","coverage_below_95pct")]
  } else {
    q<-rank_samples(expression,p$U,p$min_coverage); R<-q$ranks; meta<-q$meta
  }
  for(j in seq_len(ncol(R))) R[!is.finite(R[,j]),j]<-p$median_rank[j]
  genes<-unique(c(p$pairs$symbol1,p$pairs$symbol2))
  Z<-sweep(sweep(R[,genes,drop=FALSE],2,p$mean[genes],"-"),2,p$sd[genes],"/")
  F<-(Z[,p$pairs$symbol1,drop=FALSE]+Z[,p$pairs$symbol2,drop=FALSE])/2
  if(weighted) F<-sweep(F,2,p$pairs$W,"*")
  colnames(F)<-p$feature_order
  F[meta$status!="accepted",]<-NA_real_
  list(x=F,meta=meta)
}
predict_sdds <- function(bundle,expression,already_ranked=FALSE) {
  stopifnot(identical(bundle$version,VERSION))
  f<-transform_pair_features(expression,bundle$preprocessor,bundle$weighted,already_ranked)
  ok<-f$meta$status=="accepted"; prob<-rep(NA_real_,nrow(f$x))
  if(any(ok)) {
    if(bundle$family=="xgboost")
      prob[ok]<-as.numeric(predict(bundle$model,xgboost::xgb.DMatrix(f$x[ok,,drop=FALSE])))
    else if(bundle$family=="elastic_net")
      prob[ok]<-as.numeric(predict(bundle$model,newx=f$x[ok,,drop=FALSE],s=bundle$lambda,type="response"))
    else stop("Unknown model family")
  }
  f$meta[,`:=`(drug=bundle$drug,version=bundle$version,probability_sensitive=prob)]
  f$meta
}
