# Prepare individual fit objects for the final model-comparison workflow.
# Uses the 500-start refit for K = 6; consumed by downstream checks and the Table S3 summary.
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
a <- commandArgs(TRUE)
if (!length(a)) stop("Provide the final model output directory.")
out <- a[1]
fits <- readRDS(file.path(out, "fits_1to8_final.rds"))
fits[[6]] <- readRDS(file.path(out, "fit6_nrep500.rds"))
for (k in 1:8) saveRDS(fits[[k]], file.path(out, sprintf("fit_k%d.rds", k)))
cat("Prepared fit_k1.rds through fit_k8.rds; K=6 uses the 500-start refit.\n")
