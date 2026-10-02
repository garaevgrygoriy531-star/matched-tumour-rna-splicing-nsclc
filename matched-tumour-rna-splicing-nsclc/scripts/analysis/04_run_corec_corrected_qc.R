#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(glmnet); library(pROC)})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("usage: Rscript 04_run_corec_corrected_qc.R <INPUT_ROOT> <OUTPUT_ROOT>")
SOURCE <- normalizePath(args[[1]],mustWork=TRUE)
AUDIT <- args[[2]]
OUT1 <- file.path(AUDIT, "01_corec_stk11_audit")
OUT2 <- file.path(AUDIT, "02_corec_pipeline_qc")
OUT3 <- file.path(AUDIT, "03_corec_final_freeze")
dir.create(OUT1, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT2, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT3, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x,p) fwrite(x,p,sep="\t",quote=FALSE,na="")

MANIFEST <- file.path(SOURCE,"02_data_acquisition/gdc_expression/GDC_EXPRESSION_SELECTED_MANIFEST.tsv")
VARIANTS <- file.path(SOURCE,"02_data_acquisition/COREC_TARGET_GENE_MC3_VARIANTS.tsv")
CN_BASE <- file.path(SOURCE,"02_data_acquisition/gdc_ascat3_cohort/COHORT_LOCAL_CN_LONG.tsv")
EVENTS <- file.path(SOURCE,"04_coreC_frozen_events.tsv")
TF <- file.path(SOURCE,"04_coreC_functional_state/COREC_TF_ACTIVITY.tsv.gz")
PR <- file.path(SOURCE,"04_coreC_functional_state/COREC_PROGENY_ACTIVITY.tsv.gz")
OLD_COUNTS <- file.path(SOURCE,"04_coreC_functional_state/COREC_BENCHMARK_COUNTS.tsv")
OLD_PROJ <- file.path(SOURCE,"04_coreC_functional_state/COREC_R2_FUNCTIONAL_STATE_PROJECTION.tsv")
STK_PRIMARY <- file.path(OUT1,"STK11_LOCAL_CN_LONG.tsv")
STK_RETRY <- file.path(OUT1,"retry_failed2/STK11_LOCAL_CN_LONG.tsv")

events <- fread(EVENTS)
if(nrow(events)!=9L || anyDuplicated(events$event_id)) stop("frozen9 identity failure")
event_cases <- unique(events$patient_barcode)
manifest <- fread(MANIFEST)
if(anyDuplicated(manifest$case_barcode)) stop("duplicate expression cases")
manifest[,histology:=fifelse(project=="TCGA-LUAD","LUAD",fifelse(project=="TCGA-LUSC","LUSC",NA_character_))]
manifest[,sample_code:=substr(sample_barcode,14,15)]
manifest[,is_primary_tumor_rna:=sample_code=="01"]
if(anyNA(manifest$histology)) stop("unsupported histology")

tf <- fread(TF); pr <- fread(PR)
tf[,feature:=paste0("TF_",feature)]; pr[,feature:=paste0("PROGENY_",feature)]
scores <- rbind(tf[,.(case_barcode,feature,score=as.numeric(score))],pr[,.(case_barcode,feature,score=as.numeric(score))])
score_wide <- dcast(scores,case_barcode~feature,value.var="score",fill=NA_real_)
score_mat <- as.matrix(score_wide[,-1]); rownames(score_mat) <- score_wide$case_barcode
if(anyDuplicated(rownames(score_mat))) stop("duplicate activity cases")

vars <- fread(VARIANTS)
vars[,variant_sample_barcode:=substr(tumor_sample_barcode,1,16)]
lof_classes <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation")
noncoding_classes <- c("Silent","Intron","IGR","RNA","3'UTR","5'UTR","Flank","Targeted_Region")
vars[,is_lof:=variant_classification %chin% lof_classes]
vars[,is_nonsyn:=!(variant_classification %chin% noncoding_classes)]

cn <- fread(CN_BASE)
cn[,gene:=fifelse(grepl("^TP53_",gene_key),"TP53",gene_key)]
stk <- rbindlist(lapply(c(STK_PRIMARY,STK_RETRY),function(p) if(file.exists(p)) fread(p) else data.table()),fill=TRUE)
if(nrow(stk)) {
  stk <- unique(stk,by=c("case_barcode","sample_barcode","gene_key","file_id"))
  stk[,gene:="STK11"]
  cn <- rbind(cn[,.(case_barcode,sample_barcode,project,gene,total_copy_number=as.numeric(total_copy_number),file_id)],
              stk[,.(case_barcode,sample_barcode,project,gene,total_copy_number=as.numeric(total_copy_number),file_id)],fill=TRUE)
} else stop("corrected STK11 CN is absent")
cn <- cn[is.finite(total_copy_number)]
cn_gene <- cn[,.(total_copy_number=min(total_copy_number),cn_rows=.N,file_ids=paste(sort(unique(file_id)),collapse=";")),by=.(case_barcode,sample_barcode,gene)]

make_labels <- function(h,g) {
  m <- manifest[histology==h & case_barcode %chin% rownames(score_mat)]
  m <- m[!case_barcode %chin% event_cases]
  lab <- m[,.(case_barcode,expression_sample_barcode=sample_barcode,is_primary_tumor_rna)]
  vg <- vars[gene==g & case_barcode %chin% lab$case_barcode]
  if(nrow(vg)) vg <- merge(vg,lab[,.(case_barcode,expression_sample_barcode)],by="case_barcode",all.x=TRUE)
  exact <- if(nrow(vg)) vg[variant_sample_barcode==expression_sample_barcode,
    .(has_lof=any(is_lof),has_nonsyn=any(is_nonsyn),matched_variant_rows=.N),by=case_barcode] else data.table(case_barcode=character(),has_lof=logical(),has_nonsyn=logical(),matched_variant_rows=integer())
  unmatched <- if(nrow(vg)) vg[variant_sample_barcode!=expression_sample_barcode & is_nonsyn,
    .(has_unmatched_nonsyn=TRUE,unmatched_variant_samples=paste(sort(unique(variant_sample_barcode)),collapse=";")),by=case_barcode] else data.table(case_barcode=character(),has_unmatched_nonsyn=logical(),unmatched_variant_samples=character())
  lab <- merge(lab,exact,by="case_barcode",all.x=TRUE)
  lab <- merge(lab,unmatched,by="case_barcode",all.x=TRUE)
  lab[is.na(has_lof),has_lof:=FALSE]; lab[is.na(has_nonsyn),has_nonsyn:=FALSE]
  lab[is.na(has_unmatched_nonsyn),has_unmatched_nonsyn:=FALSE]
  cg <- cn_gene[gene==g]
  lab <- merge(lab,cg,by.x=c("case_barcode","expression_sample_barcode"),by.y=c("case_barcode","sample_barcode"),all.x=TRUE)
  lab[,cn_available:=is.finite(total_copy_number)]
  lab[,deep_deletion:=cn_available & total_copy_number<=0]
  lab[,amplified_or_complex:=cn_available & total_copy_number>4]
  lab[,benchmark_group:=fifelse(!is_primary_tumor_rna,"EXCLUDED_NONPRIMARY_RNA",
    fifelse(has_unmatched_nonsyn,"EXCLUDED_UNMATCHED_VARIANT_SAMPLE",
      fifelse(has_lof & cn_available,"CANONICAL_LOF",
        fifelse(!has_nonsyn & cn_available & !deep_deletion & !amplified_or_complex,"CLEAN_WT","EXCLUDED"))))]
  lab
}

prep <- function(train_raw,test_raw) {
  keep <- apply(train_raw,2,function(z) sum(is.finite(z))>=0.8*nrow(train_raw))
  train_raw <- train_raw[,keep,drop=FALSE]; test_raw <- test_raw[,keep,drop=FALSE]
  meds <- apply(train_raw,2,function(z) {v<-median(z,na.rm=TRUE); ifelse(is.finite(v),v,0)})
  for(j in seq_len(ncol(train_raw))) {
    train_raw[!is.finite(train_raw[,j]),j] <- meds[j]
    test_raw[!is.finite(test_raw[,j]),j] <- meds[j]
  }
  mus <- colMeans(train_raw); sds <- apply(train_raw,2,sd); sds[!is.finite(sds)|sds==0] <- 1
  list(train=sweep(sweep(train_raw,2,mus,"-"),2,sds,"/"),test=sweep(sweep(test_raw,2,mus,"-"),2,sds,"/"),keep=keep,median=meds,mean=mus,sd=sds)
}
strat_folds <- function(y,k,seed) {
  set.seed(seed); out<-integer(length(y))
  for(cl in sort(unique(y))) {ix<-which(y==cl); out[ix]<-sample(rep(seq_len(k),length.out=length(ix)))}
  out
}
auc_rank <- function(y,p) if(length(unique(y))<2) NA_real_ else as.numeric(pROC::auc(pROC::roc(y,p,quiet=TRUE,direction="<")))
weights_for <- function(y,balanced) if(!balanced) rep(1,length(y)) else ifelse(y==1,length(y)/(2*sum(y==1)),length(y)/(2*sum(y==0)))

nested_cv <- function(Xraw,y,seed_base,balanced=FALSE) {
  nfold <- min(5L,min(table(y))); aucs<-numeric()
  for(rep in seq_len(20L)) {
    outer <- strat_folds(y,nfold,seed_base+rep)
    for(fold in seq_len(nfold)) {
      te<-which(outer==fold); tr<-which(outer!=fold)
      pp<-prep(Xraw[tr,,drop=FALSE],Xraw[te,,drop=FALSE]); ytr<-y[tr]
      inner_k<-min(5L,min(table(ytr))); inner<-strat_folds(ytr,inner_k,seed_base+10000L+rep*100L+fold)
      w<-weights_for(ytr,balanced)
      fit<-tryCatch(cv.glmnet(pp$train,ytr,family="binomial",alpha=0,standardize=FALSE,nlambda=40L,foldid=inner,weights=w,type.measure="deviance"),error=function(e) NULL)
      if(is.null(fit) || length(unique(y[te]))<2) next
      pred<-as.numeric(predict(fit,newx=pp$test,s="lambda.1se",type="response")); a<-auc_rank(y[te],pred)
      if(is.finite(a)) aucs<-c(aucs,a)
    }
  }
  list(aucs=aucs,mean=if(length(aucs))mean(aucs) else NA_real_,median=if(length(aucs))median(aucs) else NA_real_,nfold=nfold)
}

fit_project <- function(Xraw,y,event_raw,seed_base,balanced=FALSE) {
  pp<-prep(Xraw,event_raw); inner_k<-min(5L,min(table(y))); inner<-strat_folds(y,inner_k,seed_base+50000L)
  w<-weights_for(y,balanced)
  fit<-cv.glmnet(pp$train,y,family="binomial",alpha=0,standardize=FALSE,nlambda=40L,foldid=inner,weights=w,type.measure="deviance")
  train_p<-as.numeric(predict(fit,newx=pp$train,s="lambda.1se",type="response"))
  event_p<-as.numeric(predict(fit,newx=pp$test,s="lambda.1se",type="response"))
  co<-as.numeric(coef(fit,s="lambda.1se"))[-1]
  list(train=train_p,event=event_p,lambda=fit$lambda.1se,features=colnames(pp$train),coef=co)
}
state_for <- function(p,lp,wp,benchmark_status) {
  if(benchmark_status=="BENCHMARK_WEAK") return("BENCHMARK_WEAK")
  lq<-quantile(lp,c(.10,.90),na.rm=TRUE); wq<-quantile(wp,c(.10,.90),na.rm=TRUE)
  if(p>=lq[1] && p<=lq[2] && (p<wq[1] || p>wq[2])) "CANONICAL_LOF_LIKE" else if(p>=wq[1] && p<=wq[2] && (p<lq[1] || p>lq[2])) "WT_LIKE" else "INTERMEDIATE"
}

# Run and persist the STK11 flow before any model fitting. This fail-fast gate
# prevents a presentation-layer error from wasting the nested-CV computation.
stk_pre <- make_labels("LUAD","STK11")
all_luad<-manifest[histology=="LUAD"]
flow<-data.table(
  stage=c("ALL_LUAD_EXPRESSION","EXCLUDE_FROZEN_R2_CASES","PRIMARY_TUMOR_RNA_ONLY","NO_STK11_NONSYN_MUTATION","NO_STK11_TRUNCATING_OR_SPLICE_MUTATION","EXACT_MATCHED_STK11_CN_AVAILABLE","NO_DEEP_DELETION","NO_COMPLEX_CN","NO_OTHER_EXCLUSION","FINAL_CLEAN_WT"),
  n_remaining=integer(10),rule=c("all selected LUAD expression cases","exclude all frozen9 event cases","sample code 01 only","no exact-sample STK11 nonsynonymous MC3 call","nested confirmation; truncating/splice are nonsynonymous","same case+16-character sample barcode and finite ASCAT3 CN","total CN > 0","total CN <= 4","no unmatched-sample STK11 nonsynonymous call","pre-registered clean-WT definition"))
flow$n_remaining[1]<-nrow(all_luad); flow$n_remaining[2]<-nrow(stk_pre)
s<-stk_pre[is_primary_tumor_rna==TRUE]; flow$n_remaining[3]<-nrow(s)
s<-s[has_nonsyn==FALSE]; flow$n_remaining[4]<-nrow(s); flow$n_remaining[5]<-nrow(s)
s<-s[cn_available==TRUE]; flow$n_remaining[6]<-nrow(s)
s<-s[deep_deletion==FALSE]; flow$n_remaining[7]<-nrow(s)
s<-s[amplified_or_complex==FALSE]; flow$n_remaining[8]<-nrow(s)
s<-s[has_unmatched_nonsyn==FALSE]; flow$n_remaining[9]<-nrow(s); flow$n_remaining[10]<-nrow(s)
flow[,n_removed:=shift(n_remaining,fill=n_remaining[1])-n_remaining]
write_tsv(flow,file.path(OUT1,"STK11_CLEAN_WT_FLOW.tsv"))
write_tsv(stk_pre,file.path(OUT1,"STK11_CASE_LEVEL_AUDIT.tsv"))
if(Sys.getenv("COREC_FLOW_PREFLIGHT_ONLY",unset="0")=="1") {
  print(flow)
  cat("COREC_FLOW_PREFLIGHT=PASS\n")
  quit(save="no",status=0)
}

design <- data.table(histology=c("LUAD","LUAD","LUSC","LUAD","LUAD","LUSC","LUSC","LUAD"),gene=c("ERCC1","ABCC10","CUL3","KEAP1","KRAS","RB1","TP53","STK11"),biology=c("NOT_SUITABLE","NOT_SUITABLE","ELIGIBLE","ELIGIBLE","NOT_SUITABLE","ELIGIBLE","ELIGIBLE","ELIGIBLE"))
count_rows<-list(); cv_rows<-list(); proj_rows<-list(); balance_rows<-list(); label_cache<-list()
old_counts<-fread(OLD_COUNTS); old_proj<-fread(OLD_PROJ)
for(i in seq_len(nrow(design))) {
  h<-design$histology[i]; g<-design$gene[i]; lab<-if(h=="LUAD" && g=="STK11") copy(stk_pre) else make_labels(h,g); label_cache[[paste(h,g,sep="|")]]<-lab
  lof<-lab[benchmark_group=="CANONICAL_LOF"]; wt<-lab[benchmark_group=="CLEAN_WT"]
  gate<-design$biology[i]=="ELIGIBLE" && nrow(lof)>=15L && nrow(wt)>=30L
  count_rows[[length(count_rows)+1L]]<-data.table(histology=h,gene=g,biology_eligibility=design$biology[i],canonical_lof_n=nrow(lof),clean_wt_n=nrow(wt),sample_size_status=ifelse(gate,"PASS","NOT_ESTIMABLE"))
  if(!gate) next
  train_cases<-c(lof$case_barcode,wt$case_barcode); Xraw<-score_mat[train_cases,,drop=FALSE]; y<-c(rep(1,nrow(lof)),rep(0,nrow(wt)))
  cv<-nested_cv(Xraw,y,20260827L+i*1000L,FALSE)
  bstat<-ifelse(is.finite(cv$mean)&&cv$mean>=.65,"BENCHMARK_MINIMALLY_INFORMATIVE","BENCHMARK_WEAK")
  ev<-events[histology==h & gene==g & patient_barcode %chin% rownames(score_mat)]
  full<-fit_project(Xraw,y,score_mat[ev$patient_barcode,,drop=FALSE],20260827L+i*1000L,FALSE)
  lp<-full$train[seq_len(nrow(lof))]; wp<-full$train[nrow(lof)+seq_len(nrow(wt))]
  cv_rows[[length(cv_rows)+1L]]<-data.table(histology=h,gene=g,canonical_lof_n=nrow(lof),clean_wt_n=nrow(wt),mean_cv_auc=cv$mean,median_cv_auc=cv$median,outer_folds=cv$nfold,repeats=20L,inner_tuning="STRATIFIED_INNER_CV_LAMBDA_1SE",full_model_lambda=full$lambda,benchmark_status=bstat)
  for(j in seq_len(nrow(ev))) {
    p<-full$event[j]; cls<-state_for(p,lp,wp,bstat)
    top<-order(abs(full$coef),decreasing=TRUE)[seq_len(min(5,length(full$coef)))]
    proj_rows[[length(proj_rows)+1L]]<-data.table(event_id=ev$event_id[j],histology=h,gene=g,benchmark_status=bstat,LoF_state_probability=p,LoF_reference_percentile=mean(lp<=p),WT_reference_percentile=mean(wp<=p),dominant_activity_features=paste(full$features[top],collapse=";"),state_class=cls,status="SUPPORTED",reason="sample-matched CN/mutation labels; fold-contained preprocessing; nested ridge tuning")
  }
  if(g %chin% c("KEAP1","RB1")) {
    bcv<-nested_cv(Xraw,y,30260827L+i*1000L,TRUE)
    bfull<-fit_project(Xraw,y,score_mat[ev$patient_barcode,,drop=FALSE],30260827L+i*1000L,TRUE)
    blp<-bfull$train[seq_len(nrow(lof))]; bwp<-bfull$train[nrow(lof)+seq_len(nrow(wt))]
    for(j in seq_len(nrow(ev))) balance_rows[[length(balance_rows)+1L]]<-data.table(event_id=ev$event_id[j],histology=h,gene=g,unweighted_mean_cv_auc=cv$mean,balanced_mean_cv_auc=bcv$mean,auc_difference=bcv$mean-cv$mean,unweighted_probability=full$event[j],balanced_probability=bfull$event[j],unweighted_state=state_for(full$event[j],lp,wp,bstat),balanced_state=state_for(bfull$event[j],blp,bwp,ifelse(is.finite(bcv$mean)&&bcv$mean>=.65,"BENCHMARK_MINIMALLY_INFORMATIVE","BENCHMARK_WEAK")),weight_rule="n/(2*n_class)",status="SENSITIVITY_ONLY")
  }
}
counts<-rbindlist(count_rows); cvtab<-if(length(cv_rows))rbindlist(cv_rows)else data.table(); proj<-if(length(proj_rows))rbindlist(proj_rows)else data.table()
missing<-setdiff(events$event_id,proj$event_id)
if(length(missing)) for(e in missing) {er<-events[event_id==e]; proj<-rbind(proj,data.table(event_id=e,histology=er$histology,gene=er$gene,benchmark_status="NOT_ESTIMABLE",LoF_state_probability=NA_real_,LoF_reference_percentile=NA_real_,WT_reference_percentile=NA_real_,dominant_activity_features="",state_class="NOT_ESTIMABLE",status="NOT_ESTIMABLE",reason="biology or corrected sample-size gate not met"),fill=TRUE)}
proj<-merge(events[,.(event_id,ord=.I)],proj,by="event_id",all.x=TRUE); setorder(proj,ord); proj[,ord:=NULL]

implementation<-data.table(
  check=c("FROZEN9_EXCLUDED_FROM_TRAINING","R2_NOT_USED_FOR_FEATURE_SELECTION","R2_NOT_USED_FOR_HYPERPARAMETERS","PRIMARY_TUMOR_RNA_REQUIRED","MUTATION_SAMPLE_EXACT_MATCH","CN_SAMPLE_EXACT_MATCH","CN_REQUIRED_FOR_LOF_AND_WT","ORIGINAL_GLOBAL_HISTOLOGY_SCALING","CORRECTED_FOLD_CONTAINED_SCALING","ORIGINAL_LAMBDA_POLICY","CORRECTED_LAMBDA_POLICY"),
  original_status=c("PASS","PASS","PASS","FAIL","FAIL","FAIL","FAIL","GLOBAL_UNSUPERVISED_SCALING","NOT_RUN","UNSPECIFIED_DATA_DEPENDENT_MINIMUM_LAMBDA","NOT_RUN"),
  corrected_status=c("PASS","PASS","PASS","PASS","PASS","PASS","PASS","REMOVED","PASS","REPLACED","STRATIFIED_INNER_CV_LAMBDA_1SE"),
  impact=c("none","none","none","normal RNA case excluded","unmatched target mutation sample excluded","01A/01B or normal/tumor CN mismatch excluded","KEAP1/TP53/STK11 counts corrected","preprocessing leakage risk removed","training-fold parameters only","not test-label leakage but not predeclared","no outer-test leakage"))
write_tsv(implementation,file.path(OUT2,"COREC_PIPELINE_IMPLEMENTATION_AUDIT.tsv"))

scaling<-merge(cvtab,old_counts[,.(histology,gene,original_mean_cv_auc=as.numeric(mean_cv_auc),original_benchmark_status=benchmark_status)],by=c("histology","gene"),all.x=TRUE)
scaling[,auc_difference:=mean_cv_auc-original_mean_cv_auc]
scaling[,scaling_qc_status:="FOLD_CONTAINED_PREPROCESSING_PASS"]
write_tsv(scaling,file.path(OUT2,"COREC_FOLD_SCALING_QC.tsv"))
tuning<-cvtab[,.(histology,gene,outer_folds,repeats,inner_tuning,full_model_lambda,tuning_leakage_found=FALSE,tuning_qc_status="PASS_NESTED")]
write_tsv(tuning,file.path(OUT2,"COREC_RIDGE_TUNING_QC.tsv"))
write_tsv(if(length(balance_rows))rbindlist(balance_rows)else data.table(),file.path(OUT2,"COREC_BALANCED_RIDGE_SENSITIVITY.tsv"))
write_tsv(counts,file.path(OUT2,"COREC_CORRECTED_BENCHMARK_COUNTS.tsv"))
write_tsv(cvtab,file.path(OUT2,"COREC_CORRECTED_BENCHMARK_CV.tsv"))

final<-merge(events[,.(event_id,histology,gene,variant_key,junction_key,frozen_transcript,cluster_id)],old_proj[,.(event_id,original_benchmark_status=benchmark_status,original_probability=LoF_state_probability,original_state_class=state_class)],by="event_id",all.x=TRUE)
final<-merge(final,proj[,.(event_id,corrected_benchmark_status=benchmark_status,corrected_probability=LoF_state_probability,corrected_LoF_reference_percentile=LoF_reference_percentile,corrected_WT_reference_percentile=WT_reference_percentile,corrected_state_class=state_class,corrected_status=status,correction_reason=reason)],by="event_id",all.x=TRUE)
final[,final_interpretation:=fcase(corrected_state_class=="CANONICAL_LOF_LIKE","activity state concordant with canonical LoF benchmark",corrected_state_class=="INTERMEDIATE","partial/intermediate activity-state penetrance",corrected_state_class=="WT_LIKE","activity state closer to clean-WT benchmark",corrected_state_class=="BENCHMARK_WEAK","benchmark discrimination below pre-specified AUC threshold",default="not estimable under pre-specified gates")]
final[,claim_ceiling:="observational transcriptomic functional-state concordance; not molecular causality or clinical actionability"]
write_tsv(final,file.path(OUT3,"COREC_FINAL_EVENT_INTERPRETATION.tsv"))
claim<-final[,.(event_id,gene,corrected_state_class,allowed_claim=final_interpretation,prohibited_overinterpretation="do not call causal, functional validation, therapeutic sensitivity, or clinical actionability",status=corrected_status)]
write_tsv(claim,file.path(OUT3,"COREC_FINAL_CLAIM_MATRIX.tsv"))

cat("COREC_CORRECTED_COUNTS\n"); print(counts)
cat("COREC_CORRECTED_CV\n"); print(cvtab)
cat("COREC_CORRECTED_PROJECTIONS\n"); print(final[,.(event_id,gene,original_state_class,corrected_state_class,corrected_probability)])
