# Parametric bootstrap LRT (BLRT): H0 = (k-1) classes vs H1 = k classes
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
suppressMessages(library(poLCA))
a <- commandArgs(TRUE); dir <- a[1]; k <- as.integer(a[2]); B <- as.integer(a[3])
g0 <- readRDS(file.path(dir, sprintf("fit_k%d.rds", k-1))); g1 <- readRDS(file.path(dir, sprintf("fit_k%d.rds", k)))
LRobs <- 2*(g1$llik - g0$llik); N <- g0$N
f <- as.formula(paste("cbind(", paste0("Y",1:28, collapse=","), ")~1"))
best <- function(dat, K, n, start=NULL) {
  bl <- -Inf
  if (!is.null(start)) { ft <- try(poLCA(f, dat, nclass=K, probs.start=start, maxiter=20000, graphs=FALSE, verbose=FALSE), silent=TRUE)
    if (!inherits(ft,"try-error")) bl <- max(bl, ft$llik) }
  for (r in seq_len(n)) { ft <- try(poLCA(f, dat, nclass=K, nrep=1, maxiter=20000, graphs=FALSE, verbose=FALSE),silent=TRUE)
    if (!inherits(ft,"try-error")) bl <- max(bl,ft$llik) }
  if (!is.finite(bl)) stop("All BLRT fits failed for a bootstrap sample (K=",K,").")
  bl
}
checkpoint <- file.path(dir,sprintf("blrt_k%d_checkpoint.rds",k))
LR <- if(file.exists(checkpoint)) readRDS(checkpoint) else rep(NA_real_,B)
if(length(LR)!=B) stop("BLRT checkpoint B does not match requested B; move checkpoint before restarting.")
for (b in which(!is.finite(LR))) {
  set.seed(20261015L + 1000L*k + b)
  sim <- poLCA.simdata(N=N, probs=unname(g0$probs), P=g0$P)$dat
  redraw <- 0L
  while(any(vapply(sim,function(v) length(unique(v))<3L,logical(1))) && redraw<100L) {
    sim <- poLCA.simdata(N=N, probs=unname(g0$probs), P=g0$P)$dat; redraw <- redraw+1L
  }
  if(any(vapply(sim,function(v) length(unique(v))<3L,logical(1)))) stop("Could not draw a bootstrap data set containing all three response categories after 100 redraws.")
  l0 <- if(k-1==1) best(sim,1,1) else best(sim,k-1,2,start=unname(g0$probs))
  l1 <- best(sim,k,6)
  LR[b] <- 2*(l1-l0)
  saveRDS(LR,checkpoint)
  write.csv(data.frame(k=k,completed=sum(is.finite(LR)),B=B,LRobs=LRobs,max_boot=max(LR,na.rm=TRUE)),
            file.path(dir,sprintf("blrt_k%d_progress.csv",k)),row.names=FALSE)
  if(b%%10==0) { cat(sprintf("k=%d bootstrap %d/%d LRobs=%.1f max_boot=%.1f exceed=%d\n",k,b,B,LRobs,max(LR,na.rm=TRUE),sum(LR>=LRobs,na.rm=TRUE))); flush.console() }
}
p <- (sum(LR>=LRobs,na.rm=TRUE)+1)/(sum(is.finite(LR))+1)
write.csv(data.frame(k=k,B=sum(is.finite(LR)),LRobs=LRobs,boot_mean=mean(LR),boot_max=max(LR),exceed=sum(LR>=LRobs),p=p),
          file.path(dir,sprintf("blrt_k%d.csv",k)),row.names=FALSE)
cat(sprintf("DONE k=%d LRobs=%.2f p=%.4f\n",k,LRobs,p))
