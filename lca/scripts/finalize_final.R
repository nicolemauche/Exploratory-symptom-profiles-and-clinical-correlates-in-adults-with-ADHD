# Combine model-fit indices, BLRT results, and bivariate-fit summaries.
# Writes the final model-comparison table (Supplementary Table S3) and bootstrap summary.
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
a <- commandArgs(TRUE); out <- a[1]; N <- 749
fits <- lapply(1:8,function(k) readRDS(file.path(out,sprintf("fit_k%d.rds",k))))
blrt <- do.call(rbind,lapply(2:6,function(k) read.csv(file.path(out,sprintf("blrt_k%d.csv",k)))))
bvr <- read.csv(file.path(out,"aggregate_bivariate_fit.csv"))
rows <- lapply(1:8,function(k){f<-fits[[k]]; p<-f$npar; K<-k; pc<-f$predclass; po<-f$posterior; tb<-table(factor(pc,1:K))
 av_by<-if(K==1) NA_real_ else sapply(1:K,function(j) mean(po[pc==j,j]))
 occ<-if(K==1) NA_real_ else (av_by/(1-av_by))/(f$P/(1-f$P))
 att<-f$attempts
 data.frame(k=K,LL=f$llik,npar=p,AIC=-2*f$llik+2*p,AICc=-2*f$llik+2*p+2*p*(p+1)/(N-p-1),
 BIC=-2*f$llik+log(N)*p,SABIC=-2*f$llik+log((N+2)/24)*p,
 relEntropy=if(K==1) NA else 1+sum(po*log(pmax(po,1e-300)))/(N*log(K)),
 minAvePP=if(K==1) NA_real_ else min(av_by),maxAvePP=if(K==1) NA_real_ else max(av_by),
 minOCC=if(K==1) NA_real_ else min(occ),minClass=min(tb),minPct=100*min(tb)/N,
 starts=if(is.null(att)) NA_integer_ else length(att),
 starts_at_best_LL=if(is.null(att)) NA_integer_ else sum(abs(att-f$llik)<=.01,na.rm=TRUE),
 BLRT_LR=if(K %in% blrt$k) blrt$LRobs[match(K,blrt$k)] else NA,
 BLRT_p=if(K %in% blrt$k) blrt$p[match(K,blrt$k)] else NA,
 aggregate_BVR=if(K %in% bvr$k) bvr$aggregate_Pearson[match(K,bvr$k)] else NA,
 pairs_gt_9.49=if(K %in% bvr$k) bvr$pairs_gt_9.49[match(K,bvr$k)] else NA)})
write.csv(do.call(rbind,rows),file.path(out,"S3_final_R433_poLCA1601.csv"),row.names=FALSE)
b1<-read.csv(file.path(out,"bootfit_seed3301.csv")); b2<-read.csv(file.path(out,"bootfit_seed4402.csv")); bv<-c(b1$X2,b2$X2)
obs<-bvr$aggregate_Pearson[bvr$k==6]
write.csv(data.frame(observed=obs,B=length(bv),mean=mean(bv),sd=sd(bv),q95=unname(quantile(bv,.95)),
 max=max(bv),exceedances=sum(bv>=obs),p_plus1=(sum(bv>=obs)+1)/(length(bv)+1)),
 file.path(out,"aggregate_fit_bootstrap_summary.csv"),row.names=FALSE)
write.csv(data.frame(bootstrap_X2=bv),file.path(out,"aggregate_fit_bootstrap_values.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo_final.txt"))
cat("Final tables written in ",normalizePath(out),"\n")
