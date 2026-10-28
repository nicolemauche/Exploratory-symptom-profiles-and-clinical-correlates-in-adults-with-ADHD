# Aggregate bivariate Pearson fit statistics across item pairs for each fitted model.
# Writes fit summaries used in Section 3.1 and the final model-comparison table.
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
if(!requireNamespace("readxl",quietly=TRUE)) stop("Requires package readxl")
suppressMessages(library(readxl))
a <- commandArgs(TRUE); xlsx <- a[1]; out <- a[2]
d <- as.data.frame(read_excel(xlsx, sheet="IN_EX_Final"))
items <- c(paste0("ua",1:5),paste0("hy",1:3),paste0("temp",1:3),paste0("afflab",1:4),paste0("emo",1:3),paste0("des",1:5),paste0("imp",1:5))
dom <- sub("[0-9]+$","",items); X <- d[,items]; N <- nrow(X)
rows <- list()
for (k in 4:8) {
  ft <- readRDS(file.path(out, sprintf("fit_k%d.rds",k))); total <- 0; within <- 0; between <- 0; nflag <- 0; npair <- 0
  for (i in 1:27) for (j in (i+1):28) {
    O <- as.matrix(table(factor(X[[i]],1:3),factor(X[[j]],1:3))); E <- matrix(0,3,3)
    for (cl in 1:k) E <- E + ft$P[cl]*outer(ft$probs[[i]][cl,],ft$probs[[j]][cl,])
    E <- E*N; use <- E>0; stat <- sum((O[use]-E[use])^2/E[use]); total <- total+stat; npair<-npair+1
    if (stat>9.49) nflag<-nflag+1
    if(dom[i]==dom[j]) within<-within+stat else between<-between+stat
  }
  rows[[k]] <- data.frame(k=k,aggregate_Pearson=total,within_sum=within,within_mean=within/45,
    within_pairs=45,between_sum=between,between_mean=between/333,between_pairs=333,pairs_gt_9.49=nflag)
}
write.csv(do.call(rbind,rows),file.path(out,"aggregate_bivariate_fit.csv"),row.names=FALSE)
cat("Wrote aggregate_bivariate_fit.csv\n")
