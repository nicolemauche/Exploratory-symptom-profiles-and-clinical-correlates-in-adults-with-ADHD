# Run with source(file.choose()) and select this file.
required_R <- "4.3.3"; required_poLCA <- "1.6.0.1"
if(as.character(getRversion())!=required_R) stop("Requires R ",required_R,"; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!=required_poLCA)
  stop("Requires poLCA ",required_poLCA,"; current is ",if(requireNamespace("poLCA",quietly=TRUE)) as.character(packageVersion("poLCA")) else "not installed")
if(!requireNamespace("readxl",quietly=TRUE)) stop("Install readxl first.")
# Find this runner's folder when sourced in RStudio.
ofile <- NULL
for(i in seq_len(sys.nframe())) { f <- tryCatch(sys.frame(i)$ofile,error=function(e) NULL); if(!is.null(f)) ofile <- f }
if(is.null(ofile)) stop("Could not locate the run-package folder. Open and source run_final_lca.R in RStudio.")
script_dir <- dirname(normalizePath(ofile,winslash="/"))
message("Select IN_EX_Final_WRI_neu.xlsx")
xlsx <- file.choose()
d <- as.data.frame(readxl::read_excel(xlsx,sheet="IN_EX_Final"))
items <- c(paste0("ua",1:5),paste0("hy",1:3),paste0("temp",1:3),paste0("afflab",1:4),paste0("emo",1:3),paste0("des",1:5),paste0("imp",1:5))
if(nrow(d)!=749 || sum(complete.cases(d[,items]))!=749 || !("Class_neu" %in% names(d))) stop("Input does not match the verified N=749 file.")
out <- file.path(dirname(xlsx),"lca_final_R433_poLCA1601")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
log_dir <- file.path(out,"logs")
dir.create(log_dir,showWarnings=FALSE,recursive=TRUE)
# Make the RStudio library paths available to every child Rscript process.
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
step_counter <- 0L
run_child <- function(script,args) {
  step_counter <<- step_counter + 1L
  message("Running ",basename(script)); flush.console()
  cmd <- file.path(script_dir,"scripts",script)
  stem <- sprintf("%02d_%s",step_counter,tools::file_path_sans_ext(basename(script)))
  log_file <- file.path(log_dir,paste0(stem,".log"))
  err_file <- file.path(log_dir,paste0(stem,".err.log"))
  rc <- system2(rscript,args=c(shQuote(cmd),vapply(as.character(args),shQuote,character(1))),stdout=log_file,stderr=err_file)
  writeLines(c(paste("script:",script),paste("exit_code:",rc),paste("stdout:",log_file),paste("stderr:",err_file)),
              file.path(log_dir,paste0(stem,".status.txt")))
  if(!identical(as.integer(rc),0L)) stop("Step failed: ",script," (exit code ",rc,"). Check output folder before resuming.")
}
run_child("Reanalyse_LCA_ModelCheck_N749_final.R",c(xlsx,out))
run_child("prepare_fit_files_final.R",out)
run_child("starts6.R",c(xlsx,out,"1101","2500","50","100"))
run_child("starts6.R",c(xlsx,out,"2202","2500","50","100"))
run_child("eval6.R",c(out,xlsx))
run_child("check31.R",c(out,xlsx,out))
for(k in 2:6) run_child("blrt.R",c(out,k,100))
run_child("aggbvr_final.R",c(xlsx,out))
run_child("bootfit.R",c(out,out,"3301",50))
run_child("bootfit.R",c(out,out,"4402",50))
run_child("finalize_final.R",out)
cat("All final LCA checks completed. Output folder:\n",normalizePath(out),"\n")
