# Two-stage 6-class start-value search, with checkpointing and failed-fit handling.
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
if(!requireNamespace("readxl",quietly=TRUE)) stop("Requires package readxl")
suppressMessages({library(poLCA); library(readxl)})
a<-commandArgs(TRUE); d<-as.data.frame(read_excel(a[1],sheet="IN_EX_Final")); out<-a[2]
seed<-as.integer(a[3]); n1<-as.integer(a[4]); n2<-as.integer(a[5]); it1<-as.integer(a[6])
items<-c(paste0("ua",1:5),paste0("hy",1:3),paste0("temp",1:3),paste0("afflab",1:4),paste0("emo",1:3),paste0("des",1:5),paste0("imp",1:5))
f<-as.formula(paste("cbind(",paste(items,collapse=","),")~1")); X<-d[,items]
cp<-file.path(out,sprintf("stage1_seed%d_checkpoint.rds",seed))
if(file.exists(cp)) {
  z<-readRDS(cp); if(z$n1!=n1||z$it1!=it1) stop("Stage-1 checkpoint settings differ; move checkpoint before restarting.")
  if(is.null(z$rng_state)) stop("Checkpoint lacks a saved RNG state. Move it aside and restart this seed from the beginning.")
  stage1<-z$probs; ll1<-z$ll
  assign(".Random.seed",z$rng_state,envir=.GlobalEnv)
} else { stage1<-vector("list",n1); ll1<-rep(-Inf,n1) }
if(!file.exists(cp)) set.seed(seed)
for(r in which(!is.finite(ll1))) {
  ft<-try(suppressWarnings(poLCA(f,X,nclass=6,nrep=1,maxiter=it1,graphs=FALSE,verbose=FALSE,calc.se=FALSE)),silent=TRUE)
  if(!inherits(ft,"try-error")&&is.finite(ft$llik)){stage1[[r]]<-ft$probs;ll1[r]<-ft$llik}
  if(r%%500L==0L){saveRDS(list(n1=n1,it1=it1,probs=stage1,ll=ll1,rng_state=.Random.seed),cp)
    cat(sprintf("seed %d stage1 %d/%d successful=%d best=%.4f\n",seed,r,n1,sum(is.finite(ll1)),max(ll1)));flush.console()}
}
saveRDS(list(n1=n1,it1=it1,probs=stage1,ll=ll1,rng_state=.Random.seed),cp)
valid<-which(is.finite(ll1)); if(length(valid)<n2) stop("Too few successful stage-1 starts: ",length(valid),"; need ",n2)
top<-valid[order(ll1[valid],decreasing=TRUE)[seq_len(n2)]]
res<-list(); preds<-list(); probs<-list(); fits<-list()
for(j in seq_along(top)) {
  ft<-try(poLCA(f,X,nclass=6,probs.start=stage1[[top[j]]],maxiter=100000,graphs=FALSE,verbose=FALSE,calc.se=FALSE),silent=TRUE)
  if(!inherits(ft,"try-error")&&is.finite(ft$llik)){
    i<-length(res)+1L; res[[i]]<-data.frame(seed=seed,start=top[j],ll_stage1=ll1[top[j]],ll_final=ft$llik,iter=ft$numiter)
    preds[[i]]<-ft$predclass; probs[[i]]<-ft$probs; fits[[i]]<-ft
  }
  if(j%%10L==0L&&length(res)){
    R<-do.call(rbind,res); write.csv(R,file.path(out,sprintf("final_seed%d_progress.csv",seed)),row.names=FALSE)
    saveRDS(list(R=R,predclass=preds,probs=probs),file.path(out,sprintf("final_seed%d_progress.rds",seed)))
    cat(sprintf("seed %d stage2 %d/%d successful=%d best=%.4f\n",seed,j,n2,length(res),max(R$ll_final)));flush.console()
  }
}
if(!length(res)) stop("No stage-2 fits converged.")
R<-do.call(rbind,res); write.csv(R,file.path(out,sprintf("final_seed%d.csv",seed)),row.names=FALSE)
saveRDS(list(R=R,predclass=preds,probs=probs,fits=fits),file.path(out,sprintf("final_seed%d.rds",seed)))
cat(sprintf("DONE seed %d best final %.4f; within .01 of best: %d/%d\n",seed,max(R$ll_final),sum(max(R$ll_final)-R$ll_final<=.01),nrow(R)))
