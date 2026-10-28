# Evaluate the two-stage six-class start-value search against the selected 500-start fit.
# Writes stability summaries used to document model recovery in Section 3.1.
if(as.character(getRversion())!="4.3.3") stop("Requires R 4.3.3; current is ",getRversion())
if(!requireNamespace("poLCA",quietly=TRUE)||as.character(packageVersion("poLCA"))!="1.6.0.1") stop("Requires poLCA 1.6.0.1")
if(!requireNamespace("readxl",quietly=TRUE)) stop("Requires package readxl")
suppressMessages({library(readxl);library(poLCA)})
a<-commandArgs(TRUE);dir<-a[1];d<-as.data.frame(read_excel(a[2],sheet="IN_EX_Final"));ref<-as.integer(d$Class_neu);edref<-ref%in%c(1,2,4)
ari<-function(x,y){t<-table(factor(x,1:6),factor(y,1:6));n<-sum(t);ch<-function(z)z*(z-1)/2;s<-sum(ch(t));A<-sum(ch(rowSums(t)));B<-sum(ch(colSums(t)));e<-A*B/ch(n);den<-(A+B)/2-e;if(den==0)1 else (s-e)/den}
L<-lapply(c(1101,2202),function(s)readRDS(file.path(dir,sprintf("final_seed%d.rds",s))))
R<-do.call(rbind,lapply(L,function(z)z$R));pcs<-do.call(c,lapply(L,function(z)z$predclass));fit6<-readRDS(file.path(dir,"fit_k6.rds"))
best_stage2<-max(R$ll_final);best_all<-max(best_stage2,fit6$llik);R$delta_LL<-best_all-R$ll_final
R$ARI<-sapply(pcs,function(p)ari(p,ref))
R$kappaED<-sapply(pcs,function(p){m<-apply(table(factor(p,1:6),factor(ref,1:6)),1,which.max);ed<-m[as.character(p)]%in%c(1,2,4);po<-mean(ed==edref);pe<-mean(ed)*mean(edref)+mean(!ed)*mean(!edref);if(pe==1)NA_real_ else (po-pe)/(1-pe)})
write.csv(R,file.path(dir,"stability_6class_run_details.csv"),row.names=FALSE)
within1<-R$delta_LL<1
summary<-data.frame(fit500_LL=fit6$llik,best_stage2_LL=best_stage2,best_overall_LL=best_all,
 higher_stage2_than_fit500=best_stage2>fit6$llik+.01,
 stage2_runs_within_0.01=sum(R$delta_LL<=.01),stage2_runs_within_1=sum(within1),
 ARI_within1_min=min(R$ARI[within1],na.rm=TRUE),ARI_within1_max=max(R$ARI[within1],na.rm=TRUE),
 kappa_within1_min=min(R$kappaED[within1],na.rm=TRUE),kappa_within1_max=max(R$kappaED[within1],na.rm=TRUE))
write.csv(summary,file.path(dir,"stability_6class_summary.csv"),row.names=FALSE)
cat("Stage 1: 5,000 short starts; stage 2: 100 selected starts continued.\n");print(summary,row.names=FALSE)
if(best_stage2>fit6$llik+.01) stop("Two-stage search found a higher LL than fit_k6.rds. Stop before BLRT/fit bootstrap; update the selected K=6 model first.")
