args<-commandArgs(FALSE)
script<-sub("^--file=","",args[grepl("^--file=",args)])
stopifnot(length(script)==1L)
code<-dirname(normalizePath(script))
Sys.setenv(SPAGI_CODE_DIR=code)
stage<-commandArgs(TRUE)
if(!length(stage))stop("Usage: Rscript code/run_stage.R <network|prepare|pilot|train|validate|cross_assay|atlas|spatial|clinical_reference|clinical|audit|figures|documents>")
load_stage<-function(file)source(file.path(code,file),local=.GlobalEnv)
switch(stage[1],
 network={load_stage("rebuild_network.R");rebuild_network()},
 prepare={load_stage("prepare_training.R");prepare_training()},
 pilot={load_stage("train_models.R");pilot_model();load_stage("spatial_pilot.R");spatial_pilot()},
 train={load_stage("train_models.R");train_all()},
 validate={load_stage("validate_models.R");validate_models()},
 cross_assay={load_stage("cross_assay.R");cross_assay()},
 atlas={load_stage("apply_atlas_spatial.R");apply_atlas()},
 spatial={load_stage("apply_atlas_spatial.R");apply_spatial()},
 clinical_reference={load_stage("clinical_reference.R");prepare_clinical()},
 clinical={load_stage("clinical_models.R");pcr_models();survival_models()},
 audit={load_stage("final_scientific_qa.R");final_scientific_qa()},
 figures={load_stage("figures_v3.R");make_figures()},
 documents={load_stage("documents_v3.R");make_documents()},
 stop("Unknown stage: ",stage[1]))
