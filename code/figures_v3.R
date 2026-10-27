source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
suppressPackageStartupMessages({library(ggplot2);library(ggprism);library(patchwork)})
tab_v3<-function(n)fread(file.path(OUT,"tables",n),na.strings="NA")
PAL<-c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#6A3D9A","#666666")
theme_set(theme_prism(base_size=10,base_family="Arial")+
  theme(plot.title=element_text(size=11),legend.title=element_blank(),legend.position="bottom",
        axis.text=element_text(size=9),strip.text=element_text(size=10),
        plot.tag=element_text(size=13,face="bold"),plot.margin=margin(7,8,6,5)))
GROUPS<-c("Malignant_TNBC","Malignant_HER2+","Malignant_ER+","T_cells","Myeloid","Fibroblasts","Endothelial")
glabel<-function(x)gsub("_"," ",sub("^Malignant_","",x))
export_figure<-function(p,name,w,h) {
  ggsave(file.path(OUT,"figures",paste0(name,".png")),p,width=w,height=h,dpi=320,bg="white")
  ggsave(file.path(OUT,"figures",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")
}
make_figures<-function() {
  net<-readRDS(file.path(WORK,"network_v3.rds"))
  raw<-fread(file.path(ROOT,"02_GI_networks/propR_edges_all.csv"),
             select=c("group","symbol1","symbol2","retain_propr"))
  old<-raw[retain_propr==TRUE,.(status=ifelse(symbol1==symbol2,"Self-pair removed","Distinct-gene retained"))][,.N,by=.(group=raw[retain_propr==TRUE]$group,status)]
  old[,group:=factor(group,levels=rev(GROUPS))]
  p2a<-ggplot(old,aes(group,N,fill=status))+geom_col(width=.7)+
    scale_fill_manual(values=c("Distinct-gene retained"=PAL[1],"Self-pair removed"=PAL[2]))+
    coord_flip()+scale_x_discrete(labels=glabel)+labs(x=NULL,y="Group records",title="Legacy retention audit")
  ph<-copy(net$summary);ph[,group:=factor(group,levels=rev(GROUPS))]
  p2b<-ggplot(ph,aes(group,physical_edges))+geom_col(fill=PAL[3],width=.65)+
    geom_text(aes(label=physical_edges),hjust=-.2,size=3)+coord_flip()+
    scale_x_discrete(labels=glabel)+scale_y_continuous(expand=expansion(mult=c(0,.14)))+
    labs(x=NULL,y="Physical edges",title="Clean topology")
  tp<-copy(net$tops);tp[,group:=factor(group,levels=GROUPS)]
  p2c<-ggplot(tp[topS_estimable==TRUE],aes(group,topS))+
    geom_hline(yintercept=0,color="grey75",linewidth=.35)+
    geom_boxplot(width=.55,outlier.shape=NA,color=PAL[1],fill="#DCEAF2")+
    geom_point(position=position_jitter(width=.1,height=0,seed=SEED),size=.85,alpha=.55)+
    scale_x_discrete(labels=glabel)+theme(axis.text.x=element_text(angle=35,hjust=1))+
    labs(x=NULL,y="topS",title="Weight-permutation specificity")
  ev<-copy(ph)
  p2d<-ggplot(ev,aes(group,indirect_functional_retained))+
    geom_col(fill=PAL[2],width=.65)+geom_text(aes(label=indirect_functional_retained),hjust=-.2,size=3)+
    coord_flip()+scale_x_discrete(labels=glabel)+scale_y_continuous(expand=expansion(mult=c(0,.2)))+
    labs(x=NULL,y="Flagged typed records",title="Indirect evidence flags")
  export_figure((p2a|p2b)/(p2c|p2d)+plot_annotation(tag_levels="A"),"Fig2_Clean_Networks_v3",11.5,8.5)

  M<-readRDS(file.path(WORK,"modules_v3.rds"));tri<-M$tripartite;modules<-M$modules
  trc<-tri[,.N,by=.(malignant_target,tme_source)]
  p3a<-ggplot(trc,aes(factor(malignant_target,levels=GROUPS[1:3]),N,fill=tme_source))+
    geom_col(width=.65)+scale_fill_manual(values=PAL[c(1,2,3,4)])+
    scale_x_discrete(labels=glabel)+labs(x=NULL,y="Communication records",title="Tripartite selection")
  members<-copy(modules);members[,group:=factor(group,levels=GROUPS[1:3])]
  members[,gene:=factor(gene,levels=rev(sort(unique(na.omit(gene)))))]
  p3b<-ggplot(members[status=="estimable"],aes(group,gene))+
    geom_point(shape=21,size=3.2,fill=PAL[3],color="white",stroke=.3)+
    scale_x_discrete(labels=glabel,drop=FALSE)+labs(x=NULL,y=NULL,title="Module membership")+
    theme(axis.text.y=element_text(size=8),panel.grid.major=element_line(color="grey92",linewidth=.3))
  cor<-tab_v3("atlas_module_composition_tests_v3.csv")
  cor[,celltype:=factor(celltype,levels=rev(c("Fibroblasts","Myeloid","T cells")))]
  p3c<-ggplot(cor[status=="exploratory"],aes(rho,celltype,color=subtype,shape=subtype))+
    geom_vline(xintercept=0,color="grey75",linewidth=.4)+
    geom_point(position=position_dodge(width=.5),size=3)+
    scale_color_manual(values=c("ER+"=PAL[2],"TNBC"=PAL[1]))+scale_shape_manual(values=c("ER+"=17,"TNBC"=16))+
    scale_x_continuous(limits=c(-1,1),breaks=c(-1,-.5,0,.5,1))+
    labs(x="Spearman correlation",y=NULL,title="Malignant module means")
  export_figure((p3a/p3c|p3b)+plot_annotation(tag_levels="A"),"Fig3_New_ESV_Modules_v3",11.5,7.8)

  perf<-tab_v3("drug_performance_by_repeat_v3.csv");drugs<-DRUGS<-c("Paclitaxel","Olaparib","Talazoparib","Palbociclib","Doxorubicin","Afatinib","Lapatinib","Bortezomib")
  perf[,drug:=factor(drug,levels=rev(drugs))]
  p4a<-ggplot(perf[variant=="xgboost_weighted"],aes(AUROC,drug))+
    geom_vline(xintercept=.5,color="grey70",linetype=2,linewidth=.4)+
    geom_point(position=position_jitter(width=0,height=.08,seed=SEED),color=PAL[1],alpha=.55,size=1.7)+
    stat_summary(fun=mean,geom="point",shape=18,size=3.5)+
    scale_x_continuous(limits=c(.2,.9),breaks=seq(.2,.9,.2))+
    labs(x="OOF AUROC",y=NULL,title="Grouped repeated evaluation")
  weighted<-perf[variant=="xgboost_weighted",.(drug,repeat_id,weighted=AUROC)]
  unweighted<-perf[variant=="xgboost_unweighted",.(drug,repeat_id,unweighted=AUROC)]
  ab<-merge(weighted,unweighted,by=c("drug","repeat_id"))
  p4b<-ggplot(ab,aes(weighted,unweighted))+
    geom_abline(intercept=0,slope=1,color="grey65",linetype=2,linewidth=.5)+geom_point(size=2,alpha=.8,color=PAL[1])+
    coord_equal(xlim=c(.2,.85),ylim=c(.2,.85))+
    labs(x="Weighted AUROC",y="Unweighted AUROC",title="Matched-panel ablation")+
    theme(legend.position="none")
  held<-tab_v3("cross_assay_cell_line_performance_v3.csv")[
      variant=="xgboost_weighted"&evaluation=="cross_assay_heldout_cell_line"]
  held[,drug:=factor(drug,levels=rev(drugs))]
  p4c<-ggplot(held,aes(AUROC,drug,color=dataset,shape=dataset))+
    geom_vline(xintercept=.5,color="grey70",linetype=2,linewidth=.4)+
    geom_point(position=position_dodge(width=.5),size=2.5)+
    scale_color_manual(values=PAL[1:2])+scale_shape_manual(values=c(16,17))+
    scale_x_continuous(limits=c(0,1))+labs(x="Cell-line AUROC",y=NULL,title="Held-out cross-assay transfer")
  imp<-tab_v3("atlas_sDDS_SHAP_summary_v3.csv")[drug=="Paclitaxel"][order(-mean_abs_SHAP)][1:8]
  imp[,pair_id:=factor(gsub("__","–",pair_id,fixed=TRUE),levels=rev(gsub("__","–",pair_id,fixed=TRUE)))]
  p4d<-ggplot(imp,aes(mean_abs_SHAP,pair_id))+geom_col(fill=PAL[3],width=.65)+
    labs(x="Mean absolute SHAP",y=NULL,title="Paclitaxel full-fit explanation")+
    theme(axis.text.y=element_text(size=8))
  export_figure((p4a|p4b)/(p4c|p4d)+plot_annotation(tag_levels="A"),"Fig4_New_sDDS_Evaluation_v3",12,9)

  sp<-tab_v3("spatial_section_summary_v3.csv");sp[,section:=factor(section,levels=sp$section)]
  long<-melt(sp,id.vars=c("section","tissue"),measure.vars=c("ESV_low","ESV_high","sDDS_low","sDDS_high"),
             variable.name="measure",value.name="score")
  long[,measure:=as.character(measure)]
  long[,`:=`(endpoint=ifelse(startsWith(measure,"ESV"),"ESV","sDDS"),region=ifelse(endsWith(measure,"low"),"Low","High"))]
  long[,region:=factor(region,levels=c("Low","High"))]
  pairedplot<-function(ep,title,ylabel)ggplot(long[endpoint==ep],aes(region,score,group=section,color=tissue,shape=tissue))+
    geom_line(alpha=.6,linewidth=.5)+geom_point(size=2)+scale_color_manual(values=c("Cancer"=PAL[1],"Normal"=PAL[2]))+
    scale_shape_manual(values=c("Cancer"=16,"Normal"=17))+labs(x="Malignant-signature region",y=ylabel,title=title)
  p5a<-pairedplot("ESV","Regional ESV","ESV")
  p5b<-pairedplot("sDDS","Paclitaxel sensitive-class probability","Sensitive-class probability")
  dif<-rbindlist(list(sp[,.(section,tissue,endpoint="ESV",difference=ESV_low-ESV_high)],
                       sp[,.(section,tissue,endpoint="sDDS",difference=sDDS_low-sDDS_high)]))
  dif[,endpoint:=factor(endpoint,levels=c("ESV","sDDS"),labels=c("ESV","Paclitaxel prediction"))]
  p5c<-ggplot(dif,aes(difference,section,color=tissue,shape=tissue))+
    geom_vline(xintercept=0,color="grey70",linetype=2,linewidth=.4)+geom_point(size=2)+
    facet_wrap(~endpoint,scales="free_x")+scale_color_manual(values=c("Cancer"=PAL[1],"Normal"=PAL[2]))+
    scale_shape_manual(values=c("Cancer"=16,"Normal"=17))+
    scale_x_continuous(breaks=scales::breaks_pretty(n=3))+
    labs(x="Low − high",y="Section",title="Paired regional differences")+
    guides(color="none",shape="none")+
    theme(axis.text.y=element_text(size=8),panel.spacing.x=unit(.65,"cm"))
  p5d<-ggplot(sp,aes(Moran_ESV,section,color=tissue,shape=tissue))+
    geom_vline(xintercept=0,color="grey75",linewidth=.4)+geom_point(size=2)+
    scale_color_manual(values=c("Cancer"=PAL[1],"Normal"=PAL[2]))+
    scale_shape_manual(values=c("Cancer"=16,"Normal"=17))+
    labs(x="Moran's I",y="Section",title="Array 6-NN autocorrelation")+
    guides(color="none",shape="none")+
    theme(axis.text.y=element_text(size=8))
  export_figure((p5a|p5b)/(p5c|p5d)+plot_layout(guides="collect")+plot_annotation(tag_levels="A"),"Fig5_Rebuilt_Spatial_v3",11.5,9)

  pcr<-tab_v3("pCR_OOF_predictions_v3.csv")
  mp<-pcr[,.(y=y[1],p=mean(probability_pCR)),by=.(ResearchID,model)]
  rocdata<-mp[,{
    o<-pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)
    .(FPR=1-as.numeric(o$specificities),TPR=as.numeric(o$sensitivities))
  },by=model]
  p6a<-ggplot(rocdata,aes(FPR,TPR,color=model))+geom_path(linewidth=.8)+
    geom_abline(intercept=0,slope=1,color="grey70",linetype=2,linewidth=.4)+
    scale_color_manual(values=PAL[1:3])+coord_equal()+labs(x="False-positive rate",y="True-positive rate",title="Mean OOF predictions")
  pp<-tab_v3("pCR_OOF_performance_by_repeat_v3.csv")
  p6b<-ggplot(pp,aes(model,AUROC,group=repeat_id))+geom_line(color="grey75",linewidth=.45)+
    geom_point(color=PAL[1],size=2)+stat_summary(aes(group=model),fun=mean,geom="point",shape=18,size=4)+
    labs(x=NULL,y="OOF AUROC",title="MP-inclusive clinical models")
  tr<-tab_v3("pCR_frozen_transfer_performance_v3.csv")
  p6c<-ggplot(tr,aes(model,AUROC,group=scope,color=scope,shape=scope))+geom_line(linewidth=.6)+geom_point(size=2.5)+
    scale_color_manual(values=c("all"=PAL[1],"DOP"=PAL[2]),labels=c("All records (n = 100)","Durvalumab–olaparib–paclitaxel (n = 71)"))+
    scale_shape_manual(values=c("all"=16,"DOP"=17),labels=c("All records (n = 100)","Durvalumab–olaparib–paclitaxel (n = 71)"))+
    labs(x=NULL,y="Transfer AUROC",title="Frozen within-trial transfer")
  sv<-tab_v3("survival_ESV_results_v3.csv")
  sv[,label:=ifelse(cohort=="METABRIC",paste0("METABRIC OS, ",time_months," mo"),
                    paste0("TCGA ",endpoint,ifelse(adequate_EPV,""," †")))]
  sv[,label:=factor(label,levels=rev(label))]
  p6d<-ggplot(sv,aes(HR_per_SD,label))+geom_vline(xintercept=1,color="grey65",linetype=2,linewidth=.45)+
    geom_errorbar(aes(xmin=CI_lo,xmax=CI_hi),width=.18,orientation="y",color=PAL[1])+
    geom_point(size=2.5,color=PAL[1])+scale_x_log10()+labs(x="ESV hazard ratio per SD",y=NULL,title="Adjusted exploratory survival")
  export_figure((p6a|p6b)/(p6c|p6d)+plot_annotation(tag_levels="A"),"Fig6_Clinical_Rebuild_v3",11.5,8.8)
  for(ext in c("png","pdf")) {
    from<-file.path(Sys.getenv("SPAGI_FIG1_DIR",file.path(OUT,"figures")),paste0("Fig1_scRNA_Atlas.",ext))
    if(normalizePath(from)!=normalizePath(file.path(OUT,"figures",basename(from)))) stopifnot(system2("cp",c(shQuote(from),shQuote(file.path(OUT,"figures",basename(from)))))==0)
  }
  cat("New figures 2–6 exported as PNG/PDF; author-generated phase-1 Figure 1 carried forward unchanged.\n")
}
if(sys.nframe()==0L)make_figures()
