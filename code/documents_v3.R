source(file.path(Sys.getenv("SPAGI_CODE_DIR", "code"), "config.R"))
suppressPackageStartupMessages({library(officer);library(xml2)})
DOCWORK<-file.path(WORK,"docs");dir.create(DOCWORK,showWarnings=FALSE)
tab<-function(n)fread(file.path(OUT,"tables",n),na.strings="NA")
fmt<-function(x,d=3)formatC(x,format="f",digits=d)
newdoc<-function(tt){
 x<-read_docx()
 x<-body_set_default_section(x,prop_section(page_size=page_size(width=8.5,height=11),
    page_margins=page_mar(top=.8,bottom=.8,left=1,right=1)))
 body_add_fpar(x,fpar(ftext(tt,fp_text(font.family="Arial",font.size=20,bold=TRUE,color="#111111"))))
}
para<-function(x,t,size=11,bold=FALSE)body_add_fpar(x,
 fpar(ftext(t,fp_text(font.family="Arial",font.size=size,bold=bold,color="#2C2A26")),
 fp_p=fp_par(padding.bottom=7)))
heading<-function(x,t)body_add_fpar(x,fpar(ftext(t,
 fp_text(font.family="Arial",font.size=14,bold=TRUE,color="#111111")),
 fp_p=fp_par(keep_with_next=TRUE,padding.top=6,padding.bottom=6)))
addtab<-function(x,t){
 a<-as.data.frame(t)
 for(j in seq_along(a))if(is.numeric(a[[j]]))a[[j]]<-if(all(a[[j]]==floor(a[[j]]),na.rm=TRUE))as.character(a[[j]]) else fmt(a[[j]],4)
 esc<-function(z){z<-gsub("&","&amp;",z,fixed=TRUE);z<-gsub("<","&lt;",z,fixed=TRUE);gsub(">","&gt;",z,fixed=TRUE)}
 weights<-rep(1,ncol(a));weights[1]<-if(ncol(a)<=4)2 else 1.6
 widths<-round(weights/sum(weights)*9360);widths[1]<-widths[1]+9360-sum(widths)
 rowxml<-function(v,header=FALSE,i=0){
  cells<-vapply(seq_along(v),function(j)paste0('<w:tc><w:tcPr><w:tcW w:w="',widths[j],
   '" w:type="dxa"/><w:shd w:val="clear" w:fill="',if(header)"D4A04A" else if(i%%2)"F9F7F3" else "FFFFFF",
   '"/><w:tcMar><w:top w:w="80" w:type="dxa"/><w:bottom w:w="80" w:type="dxa"/><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tcMar></w:tcPr>',
   '<w:p><w:pPr><w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto"/></w:pPr><w:r><w:rPr>',
   '<w:rFonts w:ascii="Arial" w:hAnsi="Arial"/><w:sz w:val="18"/><w:color w:val="',if(header)"FFFFFF" else "2C2A26",
   '"/>',if(header)"<w:b/>" else "",'</w:rPr><w:t xml:space="preserve">',esc(as.character(v[j])),'</w:t></w:r></w:p></w:tc>'),character(1))
  paste0('<w:tr><w:trPr><w:cantSplit/>',if(header)'<w:tblHeader/>' else '', '</w:trPr>',paste(cells,collapse=""),'</w:tr>')
 }
 tbl<-paste0('<w:tbl xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:tblPr><w:tblW w:w="9360" w:type="dxa"/><w:jc w:val="center"/><w:tblLayout w:type="fixed"/></w:tblPr><w:tblGrid>',
 paste0('<w:gridCol w:w="',widths,'"/>',collapse=""),'</w:tblGrid>',rowxml(names(a),TRUE),
 paste(vapply(seq_len(nrow(a)),function(i)rowxml(unlist(a[i,],use.names=FALSE),FALSE,i),character(1)),collapse=""),'</w:tbl>')
 body_add_par(body_add_xml(x,tbl),"",style="Normal")
}
publish_doc<-function(x,name){
 local<-file.path(DOCWORK,name);print(x,target=local)
 stopifnot(system2("cp",c(shQuote(local),shQuote(file.path(OUT,name))))==0)
 z<-read_docx(file.path(OUT,name));stopifnot(nrow(docx_summary(z))>5)
 invisible(local)
}
make_documents<-function() {
 perf<-tab("drug_performance_summary_v3.csv")[variant=="xgboost_weighted"]
 pcr<-tab("pCR_OOF_summary_v3.csv")
 transfer<-tab("pCR_frozen_transfer_performance_v3.csv")
 sp<-tab("spatial_paired_statistics_v3.csv")
 surv<-tab("survival_ESV_results_v3.csv")
 modules<-tab("ESV_modules_v3.csv")
 title<-"Reconstruction and exploratory evaluation of rank-based tumour microenvironment-associated gene-pair scores in breast cancer"
 abstract<-paste(
  "Expression-derived gene-pair scores require consistent preprocessing and evaluation at the biological sampling unit.",
  "We reconstructed an author-developed breast-cancer gene-pair workflow using an existing breast cancer atlas, curated synthetic-lethality/viability annotations and frozen ligand–receptor outputs.",
  "Removing 216 self-pairs from 457 previously retained group records and repairing group-aware joins yielded 241 typed records representing 235 within-group physical edges; none passed Benjamini–Hochberg correction.",
  "Rebuilt subtype modules contained 19, 10 and one genes for TNBC, HER2-positive and ER-positive malignant cells, respectively.",
  "A common 13,357-gene single-sample rank interface supported 152 pair features and eight drug models.",
  "Grouped nested evaluation yielded mean repeat-specific out-of-fold AUROCs of 0.415–0.718; topology weighting produced identical XGBoost predictions to the matched unweighted control.",
  "Across 15 Visium sections, ESV and paclitaxel scores were lower in low malignant-signature regions in 11 and 13 sections, respectively; donor independence remained unresolved.",
  paste0("MP-inclusive I-SPY2 clinical models had out-of-fold AUROCs ",paste(fmt(pcr$AUROC),collapse=", "),
         " for baseline, drug-score and expanded models; frozen transfer remained within the same platform trial."),
  "No adjusted ESV survival endpoint passed multiplicity correction.",
  "The rebuild supports reproducible exploratory regional and clinical associations, not validated dual-gene mechanisms, treatment-specific benefit or clinical utility.")
 status<-"VERSION 3.0.0 | Public reproducibility release; final journal-submission approval remains with the authors."
 intro<-c(
  "Single-cell and spatial profiling provide complementary views of malignant and microenvironmental heterogeneity in breast cancer. The source atlas used here integrates these modalities, but reuse of processed data does not by itself provide independent biological replication [3].",
  "The author-developed workflow connects curated synthetic-lethality (SL) and synthetic-viability (SV) annotations with expression associations and tumor-microenvironment (TME) context. In this study, these labels describe prior candidate annotations; coexpression, proportionality or single-gene dependency associations do not establish a dual-perturbation phenotype. The operational extrinsic synthetic viability (ESV) score should not be interpreted as an experimentally demonstrated rescue effect.",
  "The present version addresses specific upstream failures: self-pairs, group-identity loss in joins, duplicated topology annotations and inconsistent training/application features. It also replaces cohort-dependent input scaling with a single-sample rank interface. Consequently, differences from earlier versions cannot be attributed solely to self-pair deletion.",
  "We ask whether the repaired workflow yields reproducible scores, whether drug classifiers generalize to held-out cell lines across assays, whether spatial regional differences persist under frozen region definitions, and whether derived scores add exploratory clinical information. Negative results and failed eligibility gates are retained rather than optimized away.")
 methods<-list(
  "Study inputs and fixed scope"=c(
   "We reused the existing processed breast cancer H5AD containing 94,607 cells, 18,517 annotated malignant cells and 26 atlas sample/patient labels. Figure 1 was generated during the authors' earlier phase-1 processing of the publicly available GSE176078 count matrices using the archived Scanpy, Harmony and infercnvpy workflow, and it was carried forward unchanged into the present reconstruction. Alignment, integration, cell annotation and malignancy calling were not rerun. No panel reproduces or adapts a previously published figure. The atlas publication provides source-study context rather than the denominator for this processed subset [3].",
   "Other inputs were existing GDSC1/2 and PRISM19Q4 breast cancer response records, CCLE expression, SCP542 and GSE173634 cell-line expression, 15 GSE213688 Visium sections, two I-SPY2 expression datasets, and the existing METABRIC and TCGA clinical-expression cohorts. The two I-SPY2 sources are Agilent microarrays and derive from the same platform trial; their tissue handling differs [1, 2]. GSE213688 provides the spatial source records [5]. No new cohort was introduced. Same-cohort cBioPortal retrieval was limited to missing module/reference genes and sample-type manifests. Input identifiers, exclusions and hashes are provided with the analysis.")
  ,"Network repair and topology"=c(
   "The saved propR association table was the authoritative source. We excluded empty symbols and self-pairs, canonicalized unordered pairs and required consistency of reciprocal duplicates. Candidate annotations, associations and indirect evidence were joined on group, canonical symbols and pair type. Joins preserved row keys; group identity was never discarded. SL/SV dual annotations remained in typed records but each group-specific physical pair had one topology edge, weighted by its maximum retained rho.",
   "We retained the original exploratory rho > 0.3 and 200-permutation empirical P < 0.05 filter, with minimum attainable P = 1/201. We separately computed within-group and global BH sensitivity, without relaxing thresholds to restore discoveries. These cell-expression association P values are not patient-level network inference.",
   "Node strength specificity (topS) was recalculated using 1,000 fixed-topology weight permutations, seed 20261026. Undefined or numerically degenerate null SD yielded undefined topS and no hub designation. Hubs were the top decile among estimable nodes or SV-annotated nodes at or above the group median. Each cross-malignant-subtype drug-pair weight was W = max(rho) × log[1 + max|topS_a| + max|topS_b|]. Undefined node contributions were set to zero and flagged; nonpositive/nonfinite weights were ineligible.")
  ,"Rebuilt modules and atlas expression"=c(
   "We reused saved CellChat communication outputs at P < 0.01 and KEGG 2021 membership. Ligands satisfied the original marker-or-network-node rule; receptors satisfied the network-node-or-pathway-co-membership rule against the rebuilt hubs. CellChat infers communication from expression and prior ligand–receptor information; the resulting records are not experimentally verified signaling edges [4]. Receptor–hub pathway co-membership was not treated as directed causality.",
   "Subtype ESV modules were rebuilt without restoring old genes. Atlas raw values were existing log1p(CP10k); linear values were recovered with expm1. The standardized, negative-valued X matrix was not used to build the new reference. Atlas ESV used atlas-wide per-gene scaling of raw linear expression and patient-label ligand weights, then was summarized only among corresponding subtype malignant cells. For composition correlations, unweighted malignant module means—not TME-weighted ESV—were related to broad fibroblast, myeloid and T-cell fractions. Nine subtype-by-lineage tests were prespecified; n < 5 was not tested and BH covered the six estimable tests.")
  ,"Common drug-score interface"=c(
   "GENCODE v44 mapped version-stripped Ensembl identifiers to symbols, with ambiguities excluded and audited. Duplicate symbols were summed for linear values and median-collapsed for log/microarray values. A 13,357-gene universe U was fixed from cross-assay measurability only, without response or target-distribution selection. All 152 eligible positive-weight distinct-gene pairs in U were used.",
   "Within each sample, observed genes in U were ranked as (average tied rank − 0.5)/number observed. Measured zero expression was not missing. Coverage below 95% rejected prediction; remaining missing ranks were imputed with training-fold gene medians. Gene means/SDs and zero-variance exclusions were learned in the training fold. Pair features were W(z_a + z_b)/2. The same fit_preprocessor, transform_pair_features and predict_sdds implementation was used for training, validation and all applications. Target samples never changed U, weights, centering, scaling or feature order.",
   "Models output the probability of the relative assay-defined sensitive class, not a drug effect or patient-benefit probability. Single-sample/batch, input-order, saved-bundle reload and other-target-sample perturbation invariance were tested.")
  ,"Drug labels, nested evaluation and cross-assay transfer"=c(
   "Only Model metadata breast cell lines with available expression entered training (67 lines before drug-specific eligibility). For each dataset and drug, tertiles were frozen using all identifiable breast response reference records: values strictly below the lower tertile were sensitive, strictly above the upper tertile resistant, and equality/middle values intermediate. PRISM records were averaged over available doses/screens before labeling. Conflicting non-intermediate labels excluded the ModelID–drug pair. Labels are relative reference definitions, not clinical thresholds.",
   "Eight prespecified drugs were evaluated by ModelID-isolated stratified outer five-fold cross-validation repeated five times, with inner three-fold tuning. All preprocessing was refit inside each inner/outer training split. Elastic-net logistic regression used alpha 0.5, training-fold glmnet standardization and 50 training-fold-local lambda fractions from 1 to 10^-4 of the path maximum; a one-standard-error rule favored stronger regularization [7].",
   "Primary XGBoost models used binary logistic loss, eta 0.05, subsample 0.8, colsample_bytree 0.8, min_child_weight 3 and lambda 1, with max_depth 1/2/3 × 100/200 rounds. Inner mean log-loss selected the model; ties favored shallower/fewer-round fits. No test-set early stopping or synthetic sample augmentation was used. Binary logistic predictions have the probability interpretation specified by the software, conditional on the fitted label definition [6]. Final models used the same inner selection on all eligible training lines, not the best outer fold.",
   "Matched unweighted controls used identical columns, folds and primary-selected hyperparameters without retuning. Primary performance was the mean of five repeat-specific pooled OOF AUROCs. AUPRC, prevalence, Brier score, log-loss and calibration were also exported. The 2,000 ModelID bootstraps used mean OOF predictions: their conditional intervals describe that distinct estimand, not the mean repeat AUROC or full model-refitting uncertainty.",
   "SCP542 was aggregated from continuous linear expression rather than treated as integer UMI counts; GSE173634 used count pseudobulks. The unit of evaluation was the cell line. Familiar lines were predicted only by outer models that excluded that ModelID; all cells from the same line stayed grouped. Unseen eligible lines were reported separately using final models. These are held-out-cell-line cross-assay evaluations, not independent clinical validation. Comparisons were limited to the rebuilt models and fold-prevalence baseline; old local proxies were not official-method benchmarks.")
  ,"Spatial applications and regional statistics"=c(
   "All 15 sections (26,726 spots), including GEO-Normal M1, M7, M9 and M10, were retained. We reproduced the original Scanpy normalization/control-gene scoring recipe, old per-spot ESV and all old regional ESV means before replacing modules. The reproduction gate was absolute/relative tolerance 10^-5. Frozen low/high groups were malignant-signature scores at or below Q1 and at or above Q3, not histological margins, boundaries or invasive fronts.",
   "The new module union was scored with the same Scanpy recipe and frozen spatial TME multiplier. The paclitaxel predicted sensitive-class probability used the common rank interface. Low-minus-high contrasts were weighted equally by section, using two-sided exact signed-rank tests and BH over two endpoints. There were no zero or tied absolute differences. Directed row-standardized array-coordinate six-nearest-neighbor Moran's I was descriptive without significance P values. Unresolved donor-to-section mapping prevents treating 15 sections as 15 independent patients.")
  ,"Clinical reference, pCR models and frozen transfer"=c(
   "The clinical reference was rebuilt from raw linear atlas expression with frozen cell-state labels. Each patient-label/type mean contributed equally within type. Reference profiles and each bulk sample were ranked across a fixed 813-gene U_ref; normalized nonnegative least-squares coefficients were called relative composition proxies, not validated cell fractions. The clinical ESV was training-standardized module mean rank × mean(myCAF proxy, SPP1_TAM proxy), separately versioned from atlas/spatial ESV. The T-cell SL panel was the top ten genes by clean T-cell SL strength over both endpoints with alphabetical tie-breaking.",
   "GSE194040 supplied 988 valid pCR records (319 pCR). GSE173839 supplied 100 valid records (36 pCR), including 71 durvalumab–olaparib–paclitaxel records (29 pCR). ResearchIDs did not overlap, but trial independence was not established. Source documentation verified 0 = MP1 and 1 = MP2 in both datasets, so all final models included MP [1, 2].",
   "M0 comprised HR/HER2/MP; M1 added frozen CCLE-trained paclitaxel and olaparib predicted sensitive-class probabilities; M2 added clinical ESV and T-cell SL. ResearchID-stratified outer five-fold × five-repeat evaluation with inner three-fold ridge tuning refit all clinical transformations and predictor-median imputations in training folds. Forty training rows lacked module mean ranks; their derived missing predictor was imputed from training data, not zero-filled. Treatment-arm indicators were an internal sensitivity analysis. Frozen full-training coefficients and transforms were applied once to the transfer dataset, with no target recalibration, tuning or cohort-dependent scaling. Transfer metrics were descriptive summaries for all records and the durvalumab–olaparib–paclitaxel subgroup, not treatment-specific benefit tests.")
  ,"Survival and multiplicity"=c(
   "Official source sample manifests verified primary specimens and one expression sample per patient. The deterministic sample-selection and complete-case flow tables were exported. Cox models included age, stage I–IV as a factor, molecular subtype and both composition proxies. METABRIC used the available CLAUDIN subtype field; TCGA used its official SUBTYPE field. Missing outcomes/covariates were not imputed. Main HRs were per within-cohort SD of clinical ESV; raw-unit results were also exported.",
   "Schoenfeld diagnostics were prespecified. Categorical PH violations triggered stratification; ESV violation triggered ESV × log(time in months), with 12-/60-month HRs and an overall ESV likelihood-ratio test. BH covered METABRIC OS, TCGA OS and TCGA DSS overall tests. An events/coefficient-df ratio < 10 precluded adequately adjusted confirmatory interpretation. Conditional time-specific CIs did not override the overall multiplicity-adjusted test.")
 )
 results<-c(
  paste0("Clean network and modules. Of 457 legacy retained typed group records, 216 self-pairs were removed. The 241 remaining typed records represented 235 within-group physical edges. No edge passed within-group or global BH sensitivity; the retained graph is exploratory. Twenty-five retained typed records carried indirect evidence flags (24 malignant single-gene DepMap associations and one myeloid pathway overlap), not dual-perturbation validation (Fig. 2)."),
  paste0("The rebuilt tripartite selection contained 95 records. TNBC, HER2+ and ER+ modules contained 19, 10 and one genes, respectively (25-gene union; Fig. 3). The ER+ module comprised GJA1. None of the six estimable malignant module–composition correlations passed BH. TNBC fibroblast correlation was rho = −0.643, nominal P = 0.086, q = 0.513; all other estimable q values were 0.987. HER2+ had four patient labels and was not tested."),
  paste0("Grouped drug prediction. Drug-specific labeled sample sizes ranged from ",min(perf$n)," to ",max(perf$n),
   ". Mean repeat-specific XGBoost OOF AUROCs ranged from ",fmt(min(perf$AUROC))," to ",fmt(max(perf$AUROC)),
   " (Table 1; Fig. 4A). Afatinib had the highest point estimate; paclitaxel and olaparib estimates were ",
   fmt(perf[drug=="Paclitaxel",AUROC])," and ",fmt(perf[drug=="Olaparib",AUROC]),
   ". These small-sample results do not establish reliable broad drug prediction. Weighted and matched unweighted XGBoost OOF probabilities were exactly identical for all eight drugs (Fig. 4B); topology-weighting benefit is not supported."),
  "Held-out cross-assay AUROCs varied substantially across drugs and datasets (Fig. 4C; full line-level n, class counts and metrics in the supplement). Atlas application generated 148,136 scores across 18,517 malignant cells and eight drugs. Full-fit SHAP explanations were descriptive; proliferation/ribosomal/metabolic features remained prominent, so a generic cell-state explanation was not excluded (Fig. 4D).",
  paste0("Frozen spatial regions. All 15 per-spot and regional legacy ESV reproduction gates passed. New ESV was lower in the low group in ",sp[endpoint=="ESV",n_lower],"/15 sections (median low-minus-high ",
   fmt(sp[endpoint=="ESV",median_difference],6),"; P = ",fmt(sp[endpoint=="ESV",p],6),"; q = ",fmt(sp[endpoint=="ESV",q_BH],6),
   "). The paclitaxel predicted sensitive-class probability was lower in ",sp[endpoint=="sDDS",n_lower],"/15 (median ",fmt(sp[endpoint=="sDDS",median_difference],6),
   "; P = ",fmt(sp[endpoint=="sDDS",p],6),"; q = ",fmt(sp[endpoint=="sDDS",q_BH],6),
   "; Fig. 5). These remain section-level exploratory patterns; neither drug response nor independent patient replication was observed."),
  paste0("MP-inclusive pCR evaluation. Mean repeat-specific OOF AUROCs for M0/M1/M2 were ",paste(fmt(pcr$AUROC),collapse="/"),
   "; Brier scores were ",paste(fmt(pcr$Brier,4),collapse="/")," (Fig. 6A,B). Adding drug scores did not improve the baseline. The mean-OOF-prediction AUROC differences were −0.0012 for M1−M0 (conditional 95% CI −0.0092 to 0.0068) and 0.0318 for M2−M1 (0.0135 to 0.0501). These intervals hold predictions fixed and are not nested-refit uncertainty."),
  paste0("Frozen M0/M1/M2 transfer AUROCs were ",paste(fmt(transfer[scope=="all"][order(model)]$AUROC),collapse="/"),
   " in all 100 records and ",paste(fmt(transfer[scope=="DOP"][order(model)]$AUROC),collapse="/"),
   " in the 71 durvalumab–olaparib–paclitaxel records (Fig. 6C). These point estimates do not establish comparative superiority or treatment-specific benefit. Both sources belong to I-SPY2 [1, 2]."),
  "Adjusted survival. METABRIC OS included 1,454 complete cases and 821 events. Stage and subtype were stratified and ESV was time-dependent: HR/SD 1.167 (95% CI 0.959–1.421) at 12 months and 1.110 (1.014–1.216) at 60 months; overall ESV P = 0.0710, BH q = 0.2131. TCGA OS included 950 cases/125 events: HR/SD 0.884 (0.680–1.148), P = 0.3596, q = 0.5395. TCGA DSS included 932 cases/67 events: HR/SD 0.948 (0.672–1.336), P = 0.7603, q = 0.7603; events/df = 9.57 precluded adequately adjusted confirmatory interpretation (Fig. 6D). No overall endpoint passed BH.",
  "Eligibility failures and withdrawn claims. GSE114461 covered 12,614/13,357 U genes (94.44%), below the fixed 95% gate; no v3 temporal prediction was released. CD36–PCYOX1L was absent from the clean retained graph, so its prior efficacy claim was not retained. Old local DDS/Beyondcell/DREEP-style comparisons and old frozen models/OOF outputs were not mixed into v3 evidence.")
 discussion<-c(
  "The primary contribution of this rebuild is a traceable correction and a consistent scoring interface, rather than a validated clinical biomarker. Removing self-pairs and preserving group identity altered topology and modules; changing to rank inputs and retraining also altered score meaning. Version differences therefore combine several deliberate repairs.",
  "The matched ablation is particularly informative: topology-weighted and unweighted XGBoost predictions were identical. Positive per-feature rescaling did not supply additional predictive information in these fitted tree models. Elastic-net scaling also largely canceled through standardization, with numerical/tie differences documented in the tables. Neither result validates topS as an incremental predictive contribution.",
  "Spatial directionality persisted for the rebuilt scores, but the comparison is conditional on expression-defined quartiles and the retained TME multiplier. It is not an anatomical invasive-front experiment and is not an observed drug-response test. Cells and spots do not create additional independent patient replicates; unresolved donor mapping limits inferential claims.",
  "The expanded clinical model had higher internal point estimates than the MP-inclusive baseline, whereas drug-score addition alone did not. Composition proxies and expression scores can capture correlated tumor biology rather than a specific interaction mechanism. The transfer source was from the same platform trial, not an independent trial; no treatment interaction was tested. The absence of multiplicity-adjusted survival support further limits a broad prognostic claim.",
  "Other limitations include exploratory cell-level network associations with coarse permutation resolution, lack of direct paired perturbation evidence, small drug-specific cell-line counts, assay-relative labels, incomplete reference signatures, unvalidated composition proxies, generic proliferation/ribosomal signal, missing predictors requiring fold-only imputation and conditional rather than full-refit bootstrap uncertainty. A validated application would require independently sampled donors, direct dual perturbations and prospectively locked external clinical evaluation. No clinical-use recommendation follows from these analyses.")
 legends<-c(
  "Figure 1. Breast cancer single-cell atlas used in the reconstruction. This figure was generated in the authors' phase-1 workflow from publicly available GSE176078 count matrices and was carried forward unchanged into the present reconstruction; no panel was reproduced or adapted from Wu et al. [3]. A, workflow-derived embedding colored by major cell type. B, molecular-subtype labels. C, workflow-derived inferred malignancy classification. D, cell-type fractions by atlas sample/patient label. The analyzed object contains 94,607 cells and 26 labels; labels are not assumed to be 26 independent patients. Color mappings differ between the embedding and composition panel and are defined by each panel's legend.",
  "Figure 2. Clean network accounting and topology. A, legacy retained typed group records, partitioned into self-pairs removed and distinct-gene records. B, clean physical edges, with SL/SV annotations collapsed within group. C, estimable topS from 1,000 fixed-topology weight permutations; points are nodes, not independent patients. D, retained typed records carrying indirect functional-evidence flags. No retained edge passed BH sensitivity; evidence flags are not paired perturbation ground truth.",
  "Figure 3. Rebuilt tripartite selection and ESV modules. A (upper left), communication-record counts by malignant target and TME source. B (lower left), correlations of unweighted malignant module means with broad-lineage fractions for TNBC (n=8 labels) and ER+ (n=10); none passed BH across six estimable tests. HER2+ (n=4) was not tested. C (right), subtype module membership (19/10/1 genes; 25-gene union). Communication/pathway co-membership is not directed causal evidence.",
  "Figure 4. Rebuilt drug-score evaluation. A, five repeat-specific pooled OOF AUROCs per drug (points) and their mean (diamonds), with ModelID-isolated nested cross-validation. B, matched weighted versus unweighted XGBoost AUROCs across drugs/repeats; points lie on the identity line and probabilities are identical. C, held-out-cell-line pseudobulk AUROCs across two expression assays, not cell-level or independent clinical validation. D, eight highest mean absolute Paclitaxel SHAP contributions across malignant atlas cells, bias excluded; full-fit descriptive explanation only. Conditional confidence intervals for mean-OOF predictions are tabulated separately, not attached to the mean-repeat estimand.",
  "Figure 5. Rebuilt spatial scores in frozen regions. A,B, low/high malignant-signature quartile means joined within each of 15 sections for ESV and the paclitaxel predicted sensitive-class probability. C, section low-minus-high differences. D, descriptive ESV Moran's I from directed row-standardized array six-nearest-neighbor weights; no significance P values. Blue circles denote Cancer and orange triangles Normal (M1/M7/M9/M10). Exact signed-rank tests give ESV P=0.047913, q=0.047913 and predicted sensitive-class probability P=0.002014, q=0.004028 across two endpoints. Sections, not verified independent patients, are the units; regions are not histological boundaries.",
  "Figure 6. MP-inclusive clinical evaluation and exploratory survival. A, ROC curves for mean OOF predictions in 988 I-SPY2 training records. B, five repeat-specific pooled OOF AUROCs, linked by repeat, with mean diamonds. M0=HR/HER2/MP; M1 adds paclitaxel and olaparib predicted sensitive-class probabilities; M2 adds clinical ESV/T-cell SL. C, frozen within-I-SPY2 transfer AUROCs in all 100 valid records and 71 durvalumab–olaparib–paclitaxel records; no target recalibration. D, ESV HRs per cohort SD with 95% CIs. METABRIC estimates are time-dependent at 12/60 months. The dagger marks TCGA DSS events/df<10. No overall ESV endpoint passed BH; the isolated 60-month CI is not an overall confirmatory result.")
 write_table(data.table(figure=1:6,legend=legends),"figure_legends_v3.csv")
 main<-newdoc(title);main<-para(main,status,9,TRUE)
 main<-para(main,"Kun Fangᵃ,ᵇ,*, Suxiao Jiangᵃ,ᵇ, Meng Luoᵃ,ᵇ, Yinling Maᵃ,ᵇ")
 main<-para(main,"ᵃ Department of Surgery, Yinchuan Maternal and Child Health Hospital, Yinchuan 750001, China.\nᵇ Department of Surgery, Yinchuan Maternal and Child Health Hospital Affiliated to Ningxia Medical University, Yinchuan 750001, China.\n*Correspondence: Kun Fang, Department of Surgery, Yinchuan Maternal and Child Health Hospital Affiliated to Ningxia Medical University, Yinchuan 750001, China. Email: fangkun@nxmu.edu.cn",10)
 main<-para(main,"Keywords: breast cancer; gene-pair scores; tumor microenvironment; spatial transcriptomics; nested cross-validation; reproducibility.",10)
 main<-heading(main,"Abstract");main<-para(main,abstract)
 main<-heading(main,"Introduction");for(t in intro)main<-para(main,t)
 main<-heading(main,"Results");for(t in results)main<-para(main,t)
 main<-body_add_break(main)
 main<-heading(main,"Table 1. Primary grouped drug-model performance")
 main<-addtab(main,perf[,.(Drug=drug,N=n,`Mean repeat AUROC`=AUROC,`Mean repeat AUPRC`=AUPRC,Brier=Brier)])
 main<-para(main,"Five repeats of pooled OOF predictions; n is the number of cell lines. These are not training-fit metrics. Conditional mean-OOF confidence intervals and all controls are in Supplementary Information and CSVs.",9)
 main<-heading(main,"Discussion");for(t in discussion)main<-para(main,t)
 main<-heading(main,"Methods");for(h in names(methods)){main<-heading(main,h);for(t in methods[[h]])main<-para(main,t)}
 main<-heading(main,"Data availability")
 main<-para(main,"The source datasets are available from their original repositories under accession codes GSE176078 and SCP1039 for the breast cancer atlas, GSE213688 for the spatial transcriptomic dataset, and GSE194040 and GSE173839 for the I-SPY2 datasets [1–3, 5]. Additional inputs and official source links are listed in the reproducibility archive. A portable manifest records relative source locators, byte sizes and SHA-256 checksums for the exact files used. Source datasets and patient-, sample-, cell- or spot-level derivatives are not redistributed in the public code package; users should retrieve eligible inputs from the original repositories under their current access and reuse terms. The versioned reproducibility archive is publicly available at https://github.com/slavemaster88-sys/breast-cancer-gene-pair-workflow and is permanently archived on Zenodo at https://doi.org/10.5281/zenodo.23684869.")
 main<-heading(main,"Code availability")
 main<-para(main,"Version 3.0.0 of the author-developed R workflow and prediction interface, archived provenance scripts, software and seed records, original figures, portable source manifests, figure-source mappings and non-record-level aggregate results is available at https://github.com/slavemaster88-sys/breast-cancer-gene-pair-workflow and archived on Zenodo at https://doi.org/10.5281/zenodo.23684869. Author-original code is licensed under the MIT License; author-original documentation, figures and eligible aggregate tables are licensed under CC BY 4.0. Frozen model objects and record-level derivatives are not redistributed because they embed source-derived information governed by dataset-specific terms; the code supports reconstruction by users who obtain eligible inputs from the original sources.")
 main<-heading(main,"Acknowledgements")
 main<-para(main,"The authors thank Ling Xu for personal support.")
 main<-heading(main,"Funding")
 main<-para(main,"This work was supported by the National Natural Science Foundation of China (82560465); Ningxia Reproductive Disease Clinical Medical Research Center Project (2023YXGCJS004); Natural Science Foundation of Ningxia Hui Autonomous Region (2022AAC03748 and 2021AAC03523); Yinchuan Science and Technology Innovation Project (2025SF38 and 2023SF25); Ningxia Young Top Talent Training Program (NRSF-2024-106); and Yinchuan Science and Technology Innovation Team Development Program (YCHT2025064).")
 main<-heading(main,"Author contributions")
 main<-para(main,"K.F.: Conceptualization, funding acquisition, project administration, supervision, and writing—review and editing. S.J.: Methodology, resources, and supervision. M.L.: Investigation and methodology. Y.M.: Investigation, visualization, and writing—original draft. All authors reviewed and approved the manuscript.")
 main<-heading(main,"Competing interests")
 main<-para(main,"The authors declare no competing interests.")
 main<-heading(main,"Ethics approval and consent to participate")
 main<-para(main,"This study was a secondary analysis of public, de-identified datasets and involved no new participant recruitment, intervention, or collection of identifiable information. Ethics approval and consent to participate were therefore not applicable to this analysis.")
 main<-heading(main,"Consent for publication")
 main<-para(main,"Not applicable. No identifiable individual-level information or newly obtained participant images are reported.")
 main<-heading(main,"Use of generative artificial intelligence")
 main<-para(main,"OpenAI Codex was used for language editing, document restructuring, reference-style conversion, and submission-file preparation. It did not generate or analyse study data. The authors reviewed all AI-assisted edits and take full responsibility for the final content.")
 main<-heading(main,"Figure legends");for(t in legends)main<-para(main,t,10)
 localmain<-publish_doc(main,"01_Manuscript_SciRep_v3.docx")
 si<-newdoc("Supplementary Information — reconstructed breast-cancer gene-pair workflow, version 3")
 si<-para(si,status,9,TRUE)
 sections<-list(
  c("S1. Clean network and module accounting","network_record_accounting.csv; network_fullkey_association_evidence.csv; network_retained_typed.csv; network_physical_edges.csv; topology_topS_v3.csv; tripartite_network_v3.csv; ESV_modules_v3.csv","241 typed records and 235 within-group physical edges. No edge passes BH; topS null degeneracy is explicitly flagged. Module memberships follow the clean graph and fixed communication/pathway rules."),
  c("S2. Mapping, availability and relative labels","GENCODE_v44_symbol_mapping.csv; assay_gene_mapping_audit.csv; sDDS_common_gene_universe.csv; sDDS_fixed_pair_panel.csv; breast_response_source_labels_v3.csv; breast_response_cutoffs_v3.csv; breast_consensus_labels_v3.csv","U=13,357 and 152 eligible pairs. Middle-tertile/equality/conflicting labels are not eligible training outcomes. All cutoff values and contributing records are retained."),
  c("S3. Drug models and honest validation","drug_training_eligibility_v3.csv; drug_OOF_all_v3.csv; drug_performance_by_repeat_v3.csv; drug_performance_summary_v3.csv; drug_paired_model_comparisons_v3.csv; drug_OOF_calibration_v3.csv; final_model_invariance_checks_v3.csv","200 outer folds across eight drugs. Each test ModelID is absent from its training set. Mean repeat AUROC is distinct from AUROC of average OOF predictions. Conditional 2,000-bootstrap CIs do not refit models."),
  c("S4. Cross-assay cell-line evaluation","cross_assay_cell_line_mapping_v3.csv; cross_assay_heldout_predictions_by_fold_v3.csv; cross_assay_cell_line_predictions_v3.csv; cross_assay_cell_line_performance_v3.csv","Pseudobulk line is the unit. Familiar lines are scored only by models excluding their ModelID; unseen line cases are separate. Single-class scopes retain Brier/log-loss but do not define AUROC/AUPRC."),
  c("S5. Atlas scores and feature contributions","atlas_ESV_cells_v3.csv; atlas_ESV_patient_subtype_v3.csv; atlas_module_composition_tests_v3.csv; atlas_sDDS_cells_v3.csv; atlas_sDDS_patient_subtype_v3.csv; atlas_sDDS_SHAP_summary_v3.csv","ESV summaries are restricted to corresponding malignant subtypes. Composition tests use unweighted module means. Full-fit SHAP is descriptive, not evidence against proliferation confounding."),
  c("S6. Spatial reproduction and paired statistics","spatial_all_reproduction_gates_v3.csv; spatial_legacy_regional_mean_reproduction_v3.csv; spatial_*_scores_and_frozen_regions_v3.csv; spatial_section_summary_v3.csv; spatial_paired_statistics_v3.csv","All old per-spot and regional ESV means reproduced. Regions are frozen low/high signature quartiles. Tests are exact with no zero/tied absolute differences and BH across two endpoints; donor independence unresolved."),
  c("S7. Clinical reference and sample provenance","clinical_rank_reference_basis_v3.csv; reference_patient_type_contributors_v3.csv; clinical_reference_gene_universe_v3.csv; clinical_new_gene_availability_v3.csv; *_sample_manifest_v3.csv; clinical_feature_missingness_v3.csv","U_ref=813; reference columns are patient-label-balanced rank profiles. NNLS coefficients are relative composition proxies, not validated fractions. Targeted missing-gene retrieval used only existing cohorts."),
  c("S8. MP-inclusive pCR and frozen transfer","MP_source_coding_verification_v3.csv; pCR_locked_model_features_v3.csv; pCR_fold_training_IDs_v3.csv; pCR_OOF_predictions_v3.csv; pCR_OOF_summary_v3.csv; pCR_frozen_transfer_predictions_v3.csv; pCR_frozen_transfer_performance_v3.csv; pCR_arm_adjusted_internal_*v3.csv; pCR_paired_conditional_comparisons_v3.csv","MP is 0=MP1,1=MP2 in both sources. M0/M1/M2 all include MP. Treatment-arm sensitivity is internal only. Target coefficients/transforms are frozen without recalibration."),
  c("S9. Survival adjustment and diagnostics","survival_*_case_flow_v3.csv; survival_*_sample_selection_v3.csv; survival_complete_case_data_v3.csv; survival_Schoenfeld_diagnostics_v3.csv; survival_adjusted_coefficients_v3.csv; survival_ESV_results_v3.csv","One verified primary sample per patient; complete cases, age/stage/subtype and two proxies. PH adaptations were prespecified. TCGA DSS has low events/df; overall BH results are negative."),
  c("S10. Failed gates, version changes and audit","temporal_coverage_gate_v3.csv; historical_CD36_PCYOX1L_status_v3.csv; final_scientific_QA_v3.csv; input_checkpoint_verification_final_v3.csv; version_changes_v3.csv; figure_claim_manifest_v3.csv","Temporal coverage fails the fixed 95% threshold. CD36–PCYOX1L is historical only. v3 changes both the network and feature definition; earlier scores, temporal findings and proxy superiority claims are not v3 evidence."))
 for(s in sections){
  si<-heading(si,s[1]);si<-para(si,s[2],9);si<-para(si,s[3])
  if(startsWith(s[1],"S1."))si<-addtab(si,tab("network_summary_v3.csv")[,.(Group=group,Typed=retained_typed,Physical=physical_edges,Nodes=nodes)])
  if(startsWith(s[1],"S3.")){
   si<-addtab(si,perf[,.(Drug=drug,N=n,`Mean OOF AUC`=AUROC_mean_OOF,`CI lower`=AUROC_conditional_lo,`CI upper`=AUROC_conditional_hi)])
   si<-para(si,"The interval is for AUROC of average OOF predictions, not the Table 1 mean-repeat AUROC.",9)
  }
  if(startsWith(s[1],"S6."))si<-addtab(si,sp[,.(Endpoint=endpoint,N=n,Lower=n_lower,`Median low-high`=median_difference,P=p,Q=q_BH)])
  if(startsWith(s[1],"S8.")){
   si<-addtab(si,pcr[,.(Model=model,`Repeat AUC`=AUROC,`Mean OOF AUC`=AUROC_mean_OOF,Brier)])
   si<-addtab(si,transfer[,.(Scope=scope,Model=model,N=n,Events=n_pCR,AUROC,Brier)])
  }
  if(startsWith(s[1],"S9."))si<-addtab(si,unique(surv[,.(Cohort=cohort,Endpoint=endpoint,N=n,Events=events,`Events/df`=events_per_df,`Overall P`=overall_ESV_p,Q=overall_ESV_q_BH)]))
 }
 si<-heading(si,"Exact subtype memberships")
 for(g in unique(modules$group))si<-para(si,paste0(g,": ",paste(modules[group==g&status=="estimable",gene],collapse=", ")),10)
 si<-body_add_break(si);si<-heading(si,"Source-document register")
 si<-para(si,"The numeric citations used in the manuscript resolve to the following verified source records. They support metadata/method context, not validation of the present predictions.",10)
 src<-data.table(citation=1:7,
  title=c("GEO GSE194040: I-SPY2 expression and outcomes","GEO GSE173839: durvalumab–olaparib–paclitaxel I-SPY2","Wu et al. A single-cell and spatially resolved atlas of human breast cancers (2021)",
  "Jin et al. Inference and analysis of cell-cell communication using CellChat (2021)","GEO GSE213688: spatial source","XGBoost parameter documentation","glmnet introduction"),
  url=c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194040","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173839",
  "https://www.nature.com/articles/s41588-021-00911-1","https://www.nature.com/articles/s41467-021-21246-9","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213688",
  "https://xgboost.readthedocs.io/en/stable/parameter.html","https://glmnet.stanford.edu/articles/glmnet.html"))
 fwrite(src,file.path(OUT,"external_source_manifest_v3.csv"))
 for(i in seq_len(nrow(src)))si<-para(si,paste0("[",src$citation[i],"] ",src$title[i],". ",src$url[i]),9)
 publish_doc(si,"04_Supplementary_Information_v3.docx")
cover<-newdoc("Cover letter — v3 submission draft")
for(t in c("[Date — enter at submission]","The Editors, Scientific Reports","Dear Editors,",
   paste0("We propose “",title,"” for consideration."),
   "This study reports a full computational repair and exploratory evaluation of TME-associated gene-pair scores in breast cancer. We removed self-pairs, repaired group-aware joins, rebuilt topology/modules and retrained eight drug classifiers through a unified single-sample rank interface and grouped nested validation.",
   "The manuscript reports both preserved regional patterns and unsupported earlier claims. Topology weighting did not improve XGBoost predictions; several drug models performed near or below chance. Spatial comparisons remain conditional section-level summaries, the expanded pCR model is exploratory within I-SPY2, and no adjusted survival endpoint passed multiplicity correction. We do not claim direct dual-gene validation, independent-trial clinical validation or clinical utility.",
   "Version 3.0.0 of the reproducibility release includes R code, original figures, portable source manifests, source and file hashes, and eligible non-record-level aggregate results. It is available at https://github.com/slavemaster88-sys/breast-cancer-gene-pair-workflow and archived at https://doi.org/10.5281/zenodo.23684869. Frozen model objects and record-level derivatives remain in the authors' controlled internal audit archive and are not redistributed because their reuse is governed by source-specific terms.",
   "[Before submission, the authors should confirm final all-author approval, originality, prior-publication and concurrent-submission status, and any related manuscripts or editor communications.]",
   "Sincerely,","Kun Fang","Department of Surgery, Yinchuan Maternal and Child Health Hospital Affiliated to Ningxia Medical University, Yinchuan 750001, China","fangkun@nxmu.edu.cn"))cover<-para(cover,t)
 publish_doc(cover,"02_Cover_Letter_v3.docx")
 rs<-newdoc("Reporting summary — reconstructed breast-cancer gene-pair workflow, version 3 draft")
 items<-c(
  "Status. Source notes for author review, not an official completed Nature Portfolio form. Computational gates have passed; final submission-specific confirmations remain open.",
  "Design. Retrospective reanalysis of existing public processed data; no new experimental allocation, patient recruitment, randomization or blinding. No clinical-use claim.",
  "Sample sizes. Atlas 94,607 cells, 18,517 malignant cells, 26 sample/patient labels. Drug training 29–38 eligible cell lines per drug. Spatial 15 sections/26,726 spots, including four Normal sections; independent donors unresolved. pCR training 988/319 events, transfer 100/36 and durvalumab–olaparib–paclitaxel subgroup 71/29. Survival METABRIC OS 1,454/821, TCGA OS 950/125, TCGA DSS 932/67 complete cases/events.",
  "Exclusions. Empty symbols/self-pairs, inconsistent duplicates (failure gate), ambiguous mappings, nonpositive weights, fold-zero-variance features, intermediate/conflicting drug labels and <95% U coverage are explicit. Invalid pCR outcomes and incomplete survival covariates are excluded as documented. Measured zeros are not missing; no score-zero substitution. GSE114461 fails coverage and has no v3 temporal prediction.",
  "Sampling unit and replication. Cell-line CV is grouped by ModelID; pCR is grouped by ResearchID. Pseudobulk line is the cross-assay unit; cells and spots are not independent patient n. Both clinical response sources belong to I-SPY2 despite no ResearchID overlap. Atlas reuse is not independent replication.",
  "Network statistics. Original exploratory rho/P filter retained, with 200 permutations and P resolution 1/201. No edge passes BH. topS uses 1,000 fixed-topology weight permutations. Indirect DepMap/pathway flags are not paired perturbation ground truth.",
  "Predictive evaluation. Outer 5-fold × 5-repeat, inner 3-fold training-only preprocessing/tuning. Frozen external expression universe is availability-only. Final full-data fit repeats inner selection. MP is documented as 0=MP1,1=MP2 in both response sources. Frozen transfer has no target recalibration or tuning. Conditional bootstrap holds mean OOF predictions fixed; full model-refit uncertainty is not estimated.",
  "Statistics and multiplicity. Nine atlas correlations prespecified, six estimable and none BH significant. Spatial two-sided exact signed-rank tests, no zero/tied absolute differences, BH over two endpoints. Moran's I descriptive only. Survival age/stage/subtype/proxy adjustment, complete cases, PH diagnostics and prespecified stratified/time-dependent adaptations; BH across three overall ESV tests. TCGA DSS events/df=9.57, no adequately adjusted confirmatory inference.",
  "Data/code. R is primary; existing Scanpy is used only for frozen spatial scoring. Versions, parameters, seeds, portable input hashes, original figures, eligible aggregate results and public-package file hashes are included in version 3.0.0 at https://github.com/slavemaster88-sys/breast-cancer-gene-pair-workflow and https://doi.org/10.5281/zenodo.23684869. Large source data, record-level derivatives and frozen model objects are not redistributed because their reuse is governed by source-specific terms; they remain in the authors' controlled internal audit archive.",
  "Author fields. Author names, affiliations, correspondence, contributions, funding, conflicts, acknowledgements, ethics/consent and AI-use disclosure have been supplied. Originality, concurrent-submission status and final all-author approval for submission remain author-confirmation items.")
 for(t in items)rs<-para(rs,t)
 publish_doc(rs,"03_Reporting_Summary_v3.docx")
 review<-read_docx(localmain);review<-cursor_end(review)
 review<-body_add_break(review);review<-heading(review,"Embedded figures — author review")
 review<-para(review,"Figure 1 was generated in the authors' archived phase-1 workflow from public GSE176078 count matrices and carried forward unchanged; it is not copied or adapted from the Wu et al. publication. Figures 2–6 were regenerated from version-3 outputs. Each figure is followed by its complete legend.")
 images<-c("Fig1_scRNA_Atlas","Fig2_Clean_Networks_v3","Fig3_New_ESV_Modules_v3",
    "Fig4_New_sDDS_Evaluation_v3","Fig5_Rebuilt_Spatial_v3","Fig6_Clinical_Rebuild_v3")
 pngratio<-function(f){con<-file(f,"rb");on.exit(close(con));seek(con,16);wh<-readBin(con,"integer",2,size=4,endian="big");wh[2]/wh[1]}
 for(i in seq_along(images)){
  if(i>1L)review<-body_add_break(review)
  review<-heading(review,paste("Figure",i))
  f<-file.path(OUT,"figures",paste0(images[i],".png"));ratio<-pngratio(f);w<-min(6.5,6.2/ratio)
  review<-body_add_fpar(review,fpar(external_img(src=f,width=w,height=w*ratio),
       fp_p=fp_par(text.align="center",keep_with_next=TRUE)))
  review<-para(review,legends[i],9)
 }
 rf<-file.path(DOCWORK,"01_Manuscript_Embedded_Review_v3.docx");print(review,target=rf)
 stage<-file.path(WORK,"review_spacing_v3");if(dir.exists(stage))unlink(stage,recursive=TRUE);dir.create(stage)
 utils::unzip(rf,exdir=stage);xml<-read_xml(file.path(stage,"word/document.xml"));ns<-xml_ns(xml)
 for(z in xml_find_all(xml,"//w:p[w:r/w:drawing]/w:pPr/w:spacing",ns))xml_set_attr(z,"w:lineRule","auto",ns=ns)
 write_xml(xml,file.path(stage,"word/document.xml"))
 zip::zipr(rf,list.files(stage,full.names=TRUE,all.files=TRUE,no..=TRUE),root=stage)
 stopifnot(system2("cp",c(shQuote(rf),shQuote(file.path(OUT,basename(rf)))))==0)
 writeLines(c(paste0("Abstract words: ",length(strsplit(abstract,"\\s+")[[1]])),
   "Five v3 DOCX documents created. Repository, licensing, Figure 1 rights and author metadata resolved; submission-specific confirmations remain."),
   file.path(WORK,"document_checks_v3.txt"))
 cat("Five synchronized v3 Word documents created.\n")
 invisible(list(title=title,abstract=abstract,legends=legends))
}
if(sys.nframe()==0L)make_documents()
