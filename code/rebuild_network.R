source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
canonical <- function(d) {
  d<-copy(d)
  d[,`:=`(symbol1=trimws(symbol1),symbol2=trimws(symbol2))]
  d[,exclusion_reason:=fcase(is.na(symbol1)|is.na(symbol2)|symbol1==""|symbol2=="",
                            "missing_symbol",symbol1==symbol2,"self_pair",default="")]
  rev<-which(d$exclusion_reason=="" & d$symbol1>d$symbol2)
  if(length(rev)) {
    a<-d$symbol1[rev]; d[rev,symbol1:=symbol2]; d[rev,symbol2:=a]
    for(cols in list(c("expr_dep_c1","expr_dep_c2"),c("expr_dep_p1","expr_dep_p2"),
                     c("mean_eff_A","mean_eff_B"))) {
      if(all(cols %in% names(d))) {
        tmp<-d[[cols[1]]][rev]
        set(d,rev,cols[1],d[[cols[2]]][rev]); set(d,rev,cols[2],tmp)
      }
    }
  }
  d
}
KEY<-c("group","symbol1","symbol2","pair_type")
consistent_unique<-function(d,cols,label) {
  bad<-d[,lapply(.SD,uniqueN),by=KEY,.SDcols=cols]
  bad<-bad[rowSums(as.matrix(bad[,..cols])>1L)>0]
  if(nrow(bad)) {
    write_table(bad,paste0(label,"_inconsistent_duplicates.csv"))
    stop(label,": inconsistent reciprocal/repeated records; review required")
  }
  unique(d[,c(KEY,cols),with=FALSE],by=KEY)
}
rebuild_network <- function(include_atlas=TRUE) {
  clock<-proc.time()[3]
  candidate<-canonical(fread(file.path(ROOT,"02_GI_networks/CGIdb_filtered_pairs.csv")))
  assoc<-canonical(fread(file.path(ROOT,"02_GI_networks/propR_edges_all.csv")))
  excluded<-rbindlist(list(candidate[exclusion_reason!="",c(KEY,"exclusion_reason"),with=FALSE][,stage:="candidate"],
                           assoc[exclusion_reason!="",c(KEY,"exclusion_reason"),with=FALSE][,stage:="tested"]))
  write_table(excluded,"network_excluded_records.csv")
  raw_counts<-rbindlist(lapply(list(candidate=candidate,tested=assoc),function(d)
    d[,.(n_input=.N,n_missing=sum(exclusion_reason=="missing_symbol"),n_self=sum(exclusion_reason=="self_pair"),
         n_distinct_typed=uniqueN(.SD[exclusion_reason==""],by=c("symbol1","symbol2","pair_type")),
         n_retained_old=if("retain_propr"%in%names(d)) sum(retain_propr) else NA_integer_),
      by=group,.SDcols=c("symbol1","symbol2","pair_type","exclusion_reason")]),idcol="stage")
  write_table(raw_counts,"network_record_accounting.csv")
  cand<-candidate[exclusion_reason=="",.(candidate_records=.N,
      IntegrateScore=if(any(is.finite(IntegrateScore))) max(IntegrateScore,na.rm=TRUE) else NA_real_,
      candidate_sources=paste(sort(unique(source)),collapse=";"),
      candidate_tissues=paste(sort(unique(tissue)),collapse=";")),by=KEY]
  a<-consistent_unique(assoc[exclusion_reason==""],
         c("rho","p_perm","pearson_clr","p_pearson","retain_propr","retain_pearson"),"association")
  a[,retain_recalculated:=rho>.3 & p_perm<.05]
  stopifnot(identical(a$retain_propr,a$retain_recalculated))
  a[,q_BH_group:=p.adjust(p_perm,"BH"),by=group]
  a[,q_BH_all:=p.adjust(p_perm,"BH")]
  evidence_cols<-c("expr_dep_c1","expr_dep_p1","expr_dep_c2","expr_dep_p2",
                   "mean_eff_A","mean_eff_B","depmap_sl","depmap_sv","depmap_pass","icraft_pass","functional_pass")
  ev<-canonical(fread(file.path(ROOT,"02_GI_networks/functional_evidence_all.csv"),
                     select=c(KEY,evidence_cols)))
  ev<-consistent_unique(ev[exclusion_reason==""],evidence_cols,"evidence")
  # Full keys include group throughout. Saved erroneous association columns are never reused.
  joined<-merge(cand,a,by=KEY,all=TRUE,sort=TRUE)
  joined<-merge(joined,ev,by=KEY,all=TRUE,sort=TRUE)
  stopifnot(!anyDuplicated(joined[,..KEY]))
  stopifnot(nrow(joined)==nrow(unique(rbindlist(list(cand[,..KEY],a[,..KEY],ev[,..KEY])))))
  joined[,join_status:=fcase(is.na(candidate_records),"no_candidate",is.na(rho),"not_tested",default="matched")]
  write_table(joined,"network_fullkey_association_evidence.csv")
  retained<-joined[retain_propr==TRUE]
  stopifnot(!any(retained$symbol1==retained$symbol2))
  write_table(retained,"network_retained_typed.csv")
  physical<-retained[,.(rho=max(rho),pair_types=paste(sort(unique(pair_type)),collapse=";"),
                         n_types=uniqueN(pair_type)),by=.(group,symbol1,symbol2)]
  write_table(physical,"network_physical_edges.csv")
  set.seed(SEED)
  tops<-rbindlist(lapply(sort(unique(physical$group)),function(g) {
    d<-physical[group==g]; genes<-sort(unique(c(d$symbol1,d$symbol2)))
    i<-match(d$symbol1,genes); j<-match(d$symbol2,genes)
    A<-sparseMatrix(i=c(i,j),j=rep(seq_len(nrow(d)),2),x=1,dims=c(length(genes),nrow(d)))
    strength<-as.numeric(A%*%d$rho)
    null<-replicate(1000L,as.numeric(A%*%sample(d$rho,length(d$rho),replace=FALSE)))
    mu<-rowMeans(null); sd0<-apply(null,1,sd)
    ok<-is.finite(sd0)&sd0>sqrt(.Machine$double.eps)*pmax(1,abs(strength))
    ts<-rep(NA_real_,length(genes)); ts[ok]<-(strength[ok]-mu[ok])/sd0[ok]
    sv<-unique(c(retained[group==g & pair_type=="SV",symbol1],retained[group==g & pair_type=="SV",symbol2]))
    q90<-if(any(ok)) quantile(ts[ok],.9,names=FALSE) else Inf
    med<-if(any(ok)) median(ts[ok]) else Inf
    data.table(group=g,gene=genes,strength=strength,perm_mean=mu,perm_sd=sd0,topS=ts,
               topS_estimable=ok,is_SV_node=genes%in%sv,
               hub=ok & ((ts>=q90)|(genes%in%sv & ts>=med)))
  }))
  write_table(tops,"topology_topS_v3.csv"); write_table(tops[hub==TRUE],"module_hubs_v3.csv")
  maxnode<-tops[startsWith(group,"Malignant_"),
     .(max_abs_topS=if(any(topS_estimable)) max(abs(topS[topS_estimable])) else 0,
       all_topS_undefined=!any(topS_estimable)),by=gene]
  pair<-retained[startsWith(group,"Malignant_"),
     .(rho_max=max(rho),groups=paste(sort(unique(group)),collapse=";"),
       pair_types=paste(sort(unique(pair_type)),collapse=";")),by=.(symbol1,symbol2)]
  for(k in 1:2) {
    z<-maxnode[match(pair[[paste0("symbol",k)]],gene)]
    pair[,(paste0("topS_",k)):=z$max_abs_topS]
    pair[,(paste0("topS_undefined_",k)):=z$all_topS_undefined]
  }
  pair[,W:=rho_max*log1p(topS_1+topS_2)]
  pair[,pair_id:=paste(symbol1,symbol2,sep="__")]
  pair[,exclusion_reason:=fifelse(!is.finite(W)|W<=0,"nonpositive_or_nonfinite_W","")]
  setorder(pair,symbol1,symbol2)
  write_table(pair,"drug_pair_weights_v3.csv")
  # T-cell panel uses BOTH endpoints of the cleaned SL graph.
  te<-retained[group=="T_cells" & pair_type=="SL",.(rho=max(rho)),by=.(symbol1,symbol2)]
  tc<-rbindlist(list(te[,.(gene=symbol1,rho)],te[,.(gene=symbol2,rho)]))[
      ,.(SL_strength=sum(rho)),by=gene][order(-SL_strength,gene)]
  tc[,selected:=seq_len(.N)<=10L]
  write_table(tc,"Tcell_SL_gene_panel.csv")
  summary<-merge(joined[,.(candidate_typed=sum(!is.na(candidate_records)),tested_typed=sum(!is.na(rho)),
      retained_typed=sum(retain_propr%in%TRUE),indirect_functional_retained=sum(retain_propr%in%TRUE & functional_pass%in%TRUE),
      BH_group_retained=sum(q_BH_group<.05,na.rm=TRUE),BH_global_retained=sum(q_BH_all<.05,na.rm=TRUE)),
      by=group],physical[,.(physical_edges=.N,nodes=uniqueN(c(symbol1,symbol2))),by=group],by="group",all=TRUE)
  write_table(summary,"network_summary_v3.csv")
  cat("Clean network summary:\n"); print(summary)
  net<-list(candidate=cand,association=a,retained=retained,physical=physical,tops=tops,pairs=pair,tcell=tc,summary=summary)
  save_checkpoint(net,"network_v3.rds")
  rebuild_modules(net,include_atlas=include_atlas)
  cat("Network and module elapsed seconds:",proc.time()[3]-clock,"\n")
  invisible(net)
}
rebuild_modules<-function(net,include_atlas=TRUE) {
  lr<-fread(file.path(ROOT,"03_TME_extrinsic_SV/cellchat_all_LR.csv"))[
      pval<.01 & source %in% c("myCAF","SPP1_TAM","CD8_Tex","Endothelial") & startsWith(target,"Malignant")]
  counts<-fread(file.path(ROOT,"03_TME_extrinsic_SV/cellchat_counts.csv.gz"))
  genes<-counts[[1]]; C<-as.matrix(counts[,-1]); rownames(C)<-genes
  meta<-fread(file.path(ROOT,"03_TME_extrinsic_SV/cellchat_meta.csv"),skip=1,header=FALSE,
               col.names=c("cell_id","cc_subset","subtype","patient","celltype"))
  stopifnot(setequal(colnames(C),meta$cell_id)); meta<-meta[match(colnames(C),cell_id)]
  markers<-lapply(setNames(c("myCAF","SPP1_TAM"),c("myCAF","SPP1_TAM")),function(s) {
    own<-meta$cc_subset==s
    delta<-rowMeans(C[,own,drop=FALSE])-rowMeans(C[,!own,drop=FALSE])
    frac<-rowMeans(C[,own,drop=FALSE]>0)
    # Original nlargest(100): stable input order for exact ties.
    eligible<-which(frac>.1 & delta>0)
    genes[head(eligible[order(-delta[eligible],eligible)],100L)]
  })
  write_table(rbindlist(lapply(names(markers),function(s)data.table(source=s,gene=markers[[s]]))),"TME_marker_genes_v3.csv")
  kegg<-jsonlite::fromJSON(Sys.getenv("SPAGI_KEGG",file.path(getwd(),"data","reference","kegg_2021.json")))
  allnodes<-unique(c(net$physical[startsWith(group,"Malignant_"),symbol1],
                     net$physical[startsWith(group,"Malignant_"),symbol2]))
  tri<-rbindlist(lapply(seq_len(nrow(lr)),function(i) {
    r<-lr[i]; hubs<-net$tops[group==r$target & hub==TRUE,gene]
    paths<-kegg[vapply(kegg,function(gs) r$receptor%in%gs,logical(1))]
    hits<-sort(intersect(hubs,unlist(paths,use.names=FALSE)))
    ismarker<-r$ligand%in%markers[[r$source]]
    isnode<-r$ligand%in%allnodes
    recnode<-r$receptor%in%c(net$physical[group==r$target,symbol1],net$physical[group==r$target,symbol2])
    if(!(ismarker|isnode) | !(recnode|length(hits)>0)) return(NULL)
    data.table(ligand=r$ligand,receptor=r$receptor,tme_source=r$source,malignant_target=r$target,
      prob=r$prob,pval=r$pval,pathway_name=r$pathway_name,
      ligand_is_tme_marker=ismarker,ligand_is_network_node=isnode,receptor_is_network_node=recnode,
      KEGG_comembership_hubs=paste(hits,collapse=";"),n_KEGG_hubs=length(hits))
  }))
  modules<-rbindlist(lapply(c("Malignant_TNBC","Malignant_HER2+","Malignant_ER+"),function(g) {
    d<-tri[malignant_target==g]
    gs<-sort(unique(c(d$receptor,unlist(strsplit(d$KEGG_comembership_hubs,";",fixed=TRUE)))))
    gs<-gs[!is.na(gs)&nzchar(gs)]
    if(length(gs)) data.table(group=g,gene=gs,status="estimable") else data.table(group=g,gene=NA_character_,status="empty_module")
  }))
  write_table(tri,"tripartite_network_v3.csv"); write_table(modules,"ESV_modules_v3.csv")
  cat("Tripartite records:",nrow(tri),"\n"); print(modules[,.(n_genes=sum(!is.na(gene)),status=unique(status)),by=group])
  obj<-list(tripartite=tri,modules=modules,markers=markers)
  if(include_atlas && nrow(tri)) {
    obs<-atlas_obs()
    need<-unique(c(modules$gene,tri$ligand)); need<-need[!is.na(need)]
    E<-atlas_subset(need,linear=TRUE)
    selected<-(obs$celltype=="Fibroblasts" & obs$tme_myCAF>=median(obs[celltype=="Fibroblasts",tme_myCAF],na.rm=TRUE)) |
      (obs$celltype=="Myeloid" & obs$tme_SPP1_TAM>=median(obs[celltype=="Myeloid",tme_SPP1_TAM],na.rm=TRUE))
    selected[is.na(selected)]<-FALSE
    lig<-intersect(unique(tri$ligand),colnames(E))
    stopifnot(length(lig)>0)
    obs[,selected_TME:=selected]; obs[,ligand_mean:=rowMeans(E[,lig,drop=FALSE])]
    weights<-obs[,.(TME_abundance=mean(selected_TME),
       TME_ligand_mean=if(any(selected_TME)) mean(ligand_mean[selected_TME]) else NA_real_),by=patient]
    weights[,ligand_weight:=TME_abundance*TME_ligand_mean]
    write_table(weights,"atlas_patient_ligand_weights_v3.csv")
    scores<-rbindlist(lapply(c("Malignant_TNBC","Malignant_HER2+","Malignant_ER+"),function(g) {
      gs<-modules[group==g & status=="estimable",gene]; gs<-intersect(gs,colnames(E))
      sub<-sub("^Malignant_","",g)
      if(!length(gs)) return(NULL)
      # Frozen atlas-wide gene scaling matches the original ESV recipe;
      # only the correct subtype's malignant cells are exported/aggregated.
      X<-E[,gs,drop=FALSE]; mu<-colMeans(X); sd0<-apply(X,2,sd); ok<-is.finite(sd0)&sd0>0
      if(!any(ok)) return(NULL)
      z<-sweep(sweep(X[,ok,drop=FALSE],2,mu[ok],"-"),2,sd0[ok],"/")
      idx<-which(obs$malignancy=="Malignant" & obs$subtype==sub)
      d<-obs[idx,.(cell_id,patient,subtype)]
      d[,module_score:=rowMeans(z[idx,,drop=FALSE])]
      d[,ligand_weight:=weights$ligand_weight[match(patient,weights$patient)]]
      d[,`:=`(ESV=module_score*ligand_weight,n_module_genes=sum(ok))]
      d
    }))
    write_table(scores,"atlas_ESV_cells_v3.csv")
    means<-scores[,.(n_cells=.N,module_mean=mean(module_score),ESV_mean=mean(ESV)),by=.(patient,subtype)]
    write_table(means,"atlas_ESV_patient_subtype_v3.csv")
    proportions<-merge(CJ(patient=unique(obs$patient),celltype=unique(obs$celltype)),
                       obs[,.N,by=.(patient,celltype)],by=c("patient","celltype"),all.x=TRUE)
    proportions[is.na(N),N:=0L]
    proportions[,proportion:=N/sum(N),by=patient]
    tests<-rbindlist(lapply(unique(means$subtype),function(s) {
      rbindlist(lapply(c("Fibroblasts","Myeloid","T cells"),function(ct) {
        d<-merge(means[subtype==s],proportions[celltype==ct,.(patient,proportion)],by="patient")
        d<-d[is.finite(module_mean)&is.finite(proportion)]
        if(nrow(d)<5 || sd(d$module_mean)==0 || sd(d$proportion)==0)
          return(data.table(subtype=s,celltype=ct,n=nrow(d),rho=NA_real_,p=NA_real_,status="not_estimable"))
        t<-cor.test(d$module_mean,d$proportion,method="spearman",exact=FALSE)
        data.table(subtype=s,celltype=ct,n=nrow(d),rho=unname(t$estimate),p=t$p.value,status="exploratory")
      }))
    }))
    tests[,q_BH:=p.adjust(p,"BH")]; write_table(tests,"atlas_module_composition_tests_v3.csv")
    obj$weights<-weights; obj$module_scaling_note<-"atlas raw linear, atlas-wide per-gene z; subtype-restricted aggregation"
  }
  save_checkpoint(obj,"modules_v3.rds")
  invisible(obj)
}
if(sys.nframe()==0L) rebuild_network(include_atlas=!("--no-atlas"%in%commandArgs(TRUE)))
