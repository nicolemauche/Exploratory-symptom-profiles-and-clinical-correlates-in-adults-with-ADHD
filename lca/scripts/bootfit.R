# Parametric bootstrap of aggregate bivariate Pearson statistic for the six-class model
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
suppressMessages(library(poLCA))
a <- commandArgs(TRUE); dir <- a[1]; out <- a[2]; seed <- as.integer(a[3]); B <- as.integer(a[4])
g <- readRDS(file.path(dir,"fit_k6.rds")); N <- g$N; K <- 6
f <- as.formula(paste("cbind(",paste0("Y",1:28,collapse=","),")~1"))
aggX2 <- function(Xd,ft){tot<-0;for(i in 1:27)for(j in (i+1):28){O<-as.matrix(table(factor(Xd[,i],1:3),factor(Xd[,j],1:3)));E<-matrix(0,3,3)
 for(cl in 1:K)E<-E+ft$P[cl]*outer(ft$probs[[i]][cl,],ft$probs[[j]][cl,]);E<-N*E;use<-E>0;tot<-tot+sum((O[use]-E[use])^2/E[use])};tot}
fitbest <- function(dat){
 best<-try(poLCA(f,dat,nclass=K,probs.start=unname(g$probs),maxiter=100000,graphs=FALSE,verbose=FALSE),silent=TRUE)
 st<-list();ll<-numeric()
 for(r in 1:20){ft<-try(suppressWarnings(poLCA(f,dat,nclass=K,nrep=1,maxiter=100,graphs=FALSE,verbose=FALSE)),silent=TRUE)
  if(!inherits(ft,"try-error")){st[[length(st)+1]]<-ft$probs;ll<-c(ll,ft$llik)}}
 if(length(ll)){for(i in order(-ll)[seq_len(min(3L,length(ll)))]){ft<-try(poLCA(f,dat,nclass=K,probs.start=st[[i]],maxiter=100000,graphs=FALSE,verbose=FALSE),silent=TRUE)
  if(!inherits(ft,"try-error")&&(inherits(best,"try-error")||ft$llik>best$llik))best<-ft}}
 if(inherits(best,"try-error"))stop("All six-class bootstrap refits failed.");best
}
checkpoint<-file.path(out,sprintf("bootfit_seed%d_checkpoint.rds",seed))
res<-if(file.exists(checkpoint))readRDS(checkpoint)else rep(NA_real_,B)
if(length(res)!=B)stop("Checkpoint B does not match requested B; move checkpoint before restarting.")
for(b in which(!is.finite(res))){
 set.seed(seed+b)
 sim<-poLCA.simdata(N=N,probs=unname(g$probs),P=g$P)$dat
 redraw<-0L
 while(any(vapply(sim,function(v)length(unique(v))<3L,logical(1)))&&redraw<100L){sim<-poLCA.simdata(N=N,probs=unname(g$probs),P=g$P)$dat;redraw<-redraw+1L}
 if(any(vapply(sim,function(v)length(unique(v))<3L,logical(1)))) stop("Could not draw a bootstrap data set containing all three response categories after 100 redraws.")
 ft<-fitbest(sim);res[b]<-aggX2(sim,ft);saveRDS(res,checkpoint)
 write.csv(data.frame(seed=seed,completed=sum(is.finite(res)),B=B,max_boot=max(res,na.rm=TRUE)),file.path(out,sprintf("bootfit_seed%d_progress.csv",seed)),row.names=FALSE)
 if(b%%10==0){cat(sprintf("seed %d bootstrap %d/%d mean=%.1f max=%.1f\n",seed,b,B,mean(res,na.rm=TRUE),max(res,na.rm=TRUE)));flush.console()}}
write.csv(data.frame(seed=seed,b=seq_len(B),X2=res),file.path(out,sprintf("bootfit_seed%d.csv",seed)),row.names=FALSE)
cat(sprintf("DONE seed %d\n",seed))
