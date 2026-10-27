# Run after document rendering and manual figure/layout review.
source(file.path(Sys.getenv("SPAGI_CODE_DIR","code"),"config.R"))
build_archive<-function(){
  check<-fread(file.path(OUT,"tables/release_packaging_checks_v3.csv"))
  stopifnot(all(check$pass))
  # The manifest is deliberately self-excluded, with its integrity protected by the ZIP checksum.
  manifestpath<-file.path(OUT,"file_manifest_sha256.csv")
  if(file.exists(manifestpath))unlink(manifestpath)
  idxfiles<-list.files(file.path(OUT,"tables"),recursive=TRUE,full.names=TRUE)
  fwrite(data.table(file=sub(paste0(OUT,"/"),"",idxfiles,fixed=TRUE),
    bytes=as.numeric(file.info(idxfiles)$size),
    role=ifelse(grepl("/pilot/",idxfiles),"pilot_not_primary","final_v3_derived_table")),
    file.path(OUT,"table_index_v3.csv"))
  qa<-list(version=VERSION,date="2026-10-27",
    status="computational_rebuild_complete_repository_rights_confirmed_submission_approval_pending",
    final_scientific_checks=all(fread(file.path(OUT,"tables/final_scientific_QA_v3.csv"))$pass),
    input_files_verified=108,spatial_legacy_per_spot_gates=15,spatial_legacy_regional_mean_gates=15,
    exported_interface_and_model_checks=nrow(check),
    rendered_review_pages=13,rendered_supplement_pages=4,
    visual_review="New figures 2-6 and embedded render inspected; figure1 unchanged; table/figure clipping repaired",
    full_independent_retraining_repeated=FALSE,
    figure1_unchanged=TRUE,old_v2_archive_not_modified=TRUE,
    review_pdf_sha256=sha(file.path(OUT,"01_Manuscript_Embedded_Review_v3.pdf")),
    author_open_items=c("originality/concurrent submission","final all-author submission approval",
      "related manuscripts/editor communications/reviewer information"))
  jsonlite::write_json(qa,file.path(OUT,"release_validation_summary_v3.json"),pretty=TRUE,auto_unbox=TRUE)
  fs<-list.files(OUT,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
  fs<-fs[!file.info(fs)$isdir]
  stopifnot(length(fs)>100,all(file.info(fs)$size>0))
  rel<-substring(fs,nchar(OUT)+2L)
  stopifnot(!anyDuplicated(rel),!any(grepl("(^|/)\\.|\\.tmp$|~$",rel)))
  m<-data.table(file=rel,bytes=as.numeric(file.info(fs)$size),sha256=vapply(fs,sha,character(1)))
  fwrite(m,manifestpath)
  localroot<-file.path(WORK,"archive_stage_v3")
  if(dir.exists(localroot))unlink(localroot,recursive=TRUE);dir.create(localroot)
  stopifnot(system2("cp",c("-R",shQuote(OUT),shQuote(localroot)))==0)
  staged<-file.path(localroot,basename(OUT))
  stagedm<-fread(file.path(staged,"file_manifest_sha256.csv"))
  stopifnot(all(vapply(file.path(staged,stagedm$file),sha,character(1))==stagedm$sha256))
  localzip<-file.path(WORK,"SpaGI_BC_Scientific_Reports_submission_v3.zip")
  if(file.exists(localzip))unlink(localzip)
  zip::zipr(localzip,files=basename(staged),root=localroot,compression_level=6)
  zipout<-file.path(dirname(OUT),basename(localzip))
  stopifnot(system2("cp",c(shQuote(localzip),shQuote(zipout)))==0)
  stopifnot(file.info(zipout)$size>0,sha(zipout)==sha(localzip))
  writeLines(paste(sha(zipout),basename(zipout)),
    file.path(dirname(OUT),"SpaGI_BC_Scientific_Reports_submission_v3.sha256"))
  cat("Built ZIP:",basename(zipout),"\nFiles including manifest:",nrow(m)+1L,
    "\nUncompressed bytes:",sum(m$bytes)+file.info(manifestpath)$size,
    "\nZIP bytes:",file.info(zipout)$size,"\nZIP SHA256:",sha(zipout),"\n")
  invisible(zipout)
}
verify_archive_independently<-function(zipout){
  # This extraction is separate from the directory used to construct the archive.
  target<-file.path(WORK,"independent_unzip_v3")
  if(dir.exists(target))unlink(target,recursive=TRUE);dir.create(target)
  stopifnot(system2("unzip",c("-tq",shQuote(zipout)))==0)
  utils::unzip(zipout,exdir=target)
  base<-file.path(target,basename(OUT));m<-fread(file.path(base,"file_manifest_sha256.csv"))
  allfiles<-list.files(base,recursive=TRUE,all.files=TRUE,no..=TRUE)
  allfiles<-allfiles[!file.info(file.path(base,allfiles))$isdir]
  stopifnot(setequal(allfiles,c(m$file,"file_manifest_sha256.csv")))
  actual<-as.numeric(file.info(file.path(base,m$file))$size)
  hashes<-vapply(file.path(base,m$file),sha,character(1))
  stopifnot(all(actual>0),all(actual==m$bytes),all(hashes==m$sha256))
  # Old v2 archive identity is a separate preservation check, not inferred from filenames.
  v2<-file.path(dirname(OUT),"SpaGI_BC_Scientific_Reports_submission_v2.zip")
  v2expected<-"f431601f45e60b67054c5ccd77217e5bd1677bc4a120d81d35ce9e542ff81486"
  stopifnot(sha(v2)==v2expected)
  for(f in list.files(file.path(base,"models"),pattern="\\.rds$",full.names=TRUE))invisible(readRDS(f))
  receipt<-list(archive=basename(zipout),zip_bytes=as.numeric(file.info(zipout)$size),sha256=sha(zipout),
    extracted_regular_files=length(allfiles),manifest_files_verified=nrow(m),
    extra_or_missing_files=0,zero_byte_files=0,all_sizes_match=TRUE,all_sha256_match=TRUE,
    zip_integrity_test_pass=TRUE,all_model_bundles_readable=TRUE,v2_archive_unchanged=TRUE,
    manifest_self_excluded=TRUE)
  out<-file.path(dirname(OUT),"SpaGI_BC_submission_v3_archive_verification.json")
  jsonlite::write_json(receipt,out,pretty=TRUE,auto_unbox=TRUE)
  print(receipt);invisible(receipt)
}
if(sys.nframe()==0L){
  args<-commandArgs(TRUE)
  if("--verify-only"%in%args)
    verify_archive_independently(file.path(dirname(OUT),"SpaGI_BC_Scientific_Reports_submission_v3.zip"))
  else build_archive()
}
