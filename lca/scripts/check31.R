# Summarize agreement of the non-best two-stage six-class runs with the reported solution
# (adjusted Rand index, kappa for the ED-related grouping) and mean bivariate Pearson
# statistics within and between WRAADDS domains (Section 3.1).
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
if(!requireNamespace("readxl",quietly=TRUE)) stop("Requires package readxl")
suppressMessages({library(readxl);library(poLCA)})
a<-commandArgs(TRUE);sd<-a[1];d<-as.data.frame(read_excel(a[2],sheet="IN_EX_Final"));out<-a[3];ref<-as.integer(d$Class_neu);edref<-ref%in%c(1,2,4)
L<-lapply(c(1101,2202),function(s)readRDS(file.path(sd,sprintf("final_seed%d.rds",s))))
R<-do.call(rbind,lapply(L,function(z)z$R));pcs<-do.call(c,lapply(L,function(z)z$predclass));ft<-readRDS(file.path(out,"fit_k6.rds"));best<-max(c(ft$llik,R$ll_final));nb<-which(best-R$ll_final>.01)
ari<-function(x,y){t<-table(factor(x,1:6),factor(y,1:6));n<-sum(t);ch<-function(z)z*(z-1)/2;s<-sum(ch(t));A<-sum(ch(rowSums(t)));B<-sum(ch(colSums(t)));e<-A*B/ch(n);den<-(A+B)/2-e;if(den==0)1 else (s-e)/den}
kp<-function(p){m<-apply(table(factor(p,1:6),factor(ref,1:6)),1,which.max);ed<-m[as.character(p)]%in%c(1,2,4);po<-mean(ed==edref);pe<-mean(ed)*mean(edref)+mean(!ed)*mean(!edref);if(pe==1)NA_real_ else (po-pe)/(1-pe)}
items<-c(paste0("ua",1:5),paste0("hy",1:3),paste0("temp",1:3),paste0("afflab",1:4),paste0("emo",1:3),paste0("des",1:5),paste0("imp",1:5));dom<-sub("[0-9]+$","",items);X<-d[,items];N<-nrow(X);w<-numeric();b<-numeric()
for(i in 1:27)for(j in (i+1):28){O<-as.matrix(table(factor(X[[i]],1:3),factor(X[[j]],1:3)));E<-matrix(0,3,3);for(cl in 1:6)E<-E+ft$P[cl]*outer(ft$probs[[i]][cl,],ft$probs[[j]][cl,]);E<-N*E;use<-E>0;x<-sum((O[use]-E[use])^2/E[use]);if(dom[i]==dom[j])w<-c(w,x)else b<-c(b,x)}
res<-data.frame(nonbest_runs=length(nb),max_delta_LL=max(best-R$ll_final),ARI_min=if(length(nb))min(sapply(pcs[nb],ari,y=ref))else NA,
 ARI_max=if(length(nb))max(sapply(pcs[nb],ari,y=ref))else NA,kappa_min=if(length(nb))min(sapply(pcs[nb],kp))else NA,
 kappa_max=if(length(nb))max(sapply(pcs[nb],kp))else NA,within_sum=sum(w),within_mean=mean(w),within_pairs=length(w),
 between_sum=sum(b),between_mean=mean(b),between_pairs=length(b))
write.csv(res,file.path(out,"stability_and_domain_residual_summary.csv"),row.names=FALSE);print(res,row.names=FALSE)
