#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(decoupleR)})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=8L) stop("expected: RNA PROTEIN MAP CROSSWALK DOROTHEA PROGENY OUTROOT B")
RNA <- args[1]; PROTEIN <- args[2]; MAP <- args[3]; CROSSWALK <- args[4]
DOROTHEA <- args[5]; PROGENY <- args[6]; OUTROOT <- args[7]; B <- as.integer(args[8])
if(!is.finite(B) || B!=10000L) stop("formal bootstrap B must equal frozen 10000")
if(as.character(packageVersion("decoupleR"))!="2.16.0") stop("decoupleR version drift")

OUT3 <- file.path(OUTROOT,"03_cptac_reference_groups")
OUT4 <- file.path(OUTROOT,"04_cptac_rna_validation")
OUT5 <- file.path(OUTROOT,"05_cptac_proteomic_validation")
OUT6 <- file.path(OUTROOT,"06_cptac_integrated_results")
for(p in c(OUT3,OUT4,OUT5,OUT6)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,p) fwrite(x,p,sep="\t",quote=FALSE,na="NA")

map <- fread(MAP)
if(anyDuplicated(map$ensembl_gene_id_base)) stop("duplicate Ensembl mapping")
emap <- setNames(map$gene_symbol,map$ensembl_gene_id_base)

read_symbol_matrix <- function(path) {
  x <- fread(path)
  if(names(x)[1]!="idx") stop("matrix first column must be idx")
  samples <- names(x)[-1]
  if(anyDuplicated(samples)) stop("duplicate matrix samples")
  x[,ensembl_gene_id_base:=sub("\\..*$","",idx)]
  x[,symbol:=unname(emap[ensembl_gene_id_base])]
  x <- x[!is.na(symbol) & nzchar(symbol)]
  for(nm in samples) set(x,j=nm,value=suppressWarnings(as.numeric(x[[nm]])))
  x <- x[,lapply(.SD,function(v) {z<-mean(v,na.rm=TRUE); if(is.nan(z)) NA_real_ else z}),by=symbol,.SDcols=samples]
  m <- as.matrix(x[,-1]); rownames(m)<-x$symbol; storage.mode(m)<-"double"
  m
}

zscore <- function(v) {
  mu<-mean(v,na.rm=TRUE); s<-sd(v,na.rm=TRUE)
  if(!is.finite(s)||s==0) return(rep(NA_real_,length(v)))
  (v-mu)/s
}
cliffs <- function(x,y) mean(outer(x,y,"\u003e"))-mean(outer(x,y,"\u003c"))
hedges_g <- function(x,y) {
  nx<-length(x); ny<-length(y); sp<-sqrt(((nx-1)*var(x)+(ny-1)*var(y))/(nx+ny-2))
  if(!is.finite(sp)||sp==0) return(NA_real_)
  d<-(mean(x)-mean(y))/sp; j<-1-3/(4*(nx+ny)-9); j*d
}
classify <- function(med,lo,hi,delta,expected_sign) {
  signed<-expected_sign*med; slo<-min(expected_sign*lo,expected_sign*hi); shi<-max(expected_sign*lo,expected_sign*hi)
  if(is.finite(slo)&&slo>0&&is.finite(delta)&&abs(delta)>=.33) return("DIRECTIONALLY_REPLICATED")
  if(is.finite(signed)&&signed>0) return("PARTIALLY_REPLICATED")
  if(is.finite(shi)&&shi<0) return("OPPOSITE_DIRECTION")
  "NO_CLEAR_REPLICATION"
}
summarize_endpoint <- function(values,groups,endpoint,modality,expected_sign,seed,role="PRIMARY") {
  keep<-is.finite(values)&groups%chin%c("CANONICAL_LOF","CLEAN_WT")
  x<-values[keep & groups=="CANONICAL_LOF"]; y<-values[keep & groups=="CLEAN_WT"]
  if(length(x)<8L||length(y)<30L) return(data.table(modality=modality,endpoint=endpoint,role=role,n_lof=length(x),n_wt=length(y),median_lof=median(x),median_wt=median(y),median_difference=NA_real_,hedges_g=NA_real_,cliffs_delta=NA_real_,bootstrap_median_difference_lo=NA_real_,bootstrap_median_difference_hi=NA_real_,bootstrap_hedges_g_lo=NA_real_,bootstrap_hedges_g_hi=NA_real_,bootstrap_cliffs_delta_lo=NA_real_,bootstrap_cliffs_delta_hi=NA_real_,wilcoxon_p=NA_real_,expected_direction=ifelse(expected_sign>0,"POSITIVE","NEGATIVE"),replication_class="NOT_ESTIMABLE"))
  med<-median(x)-median(y); hg<-hedges_g(x,y); cd<-cliffs(x,y)
  set.seed(seed)
  boot<-replicate(B,{bx<-sample(x,length(x),replace=TRUE);by<-sample(y,length(y),replace=TRUE);c(median(bx)-median(by),hedges_g(bx,by),cliffs(bx,by))})
  ci<-apply(boot,1,quantile,probs=c(.025,.975),na.rm=TRUE,names=FALSE,type=7)
  wp<-suppressWarnings(wilcox.test(x,y,exact=FALSE,correct=FALSE)$p.value)
  data.table(modality=modality,endpoint=endpoint,role=role,n_lof=length(x),n_wt=length(y),median_lof=median(x),median_wt=median(y),median_difference=med,hedges_g=hg,cliffs_delta=cd,bootstrap_median_difference_lo=ci[1,1],bootstrap_median_difference_hi=ci[2,1],bootstrap_hedges_g_lo=ci[1,2],bootstrap_hedges_g_hi=ci[2,2],bootstrap_cliffs_delta_lo=ci[1,3],bootstrap_cliffs_delta_hi=ci[2,3],wilcoxon_p=wp,expected_direction=ifelse(expected_sign>0,"POSITIVE","NEGATIVE"),replication_class=classify(med,ci[1,1],ci[2,1],cd,expected_sign))
}

cw <- fread(CROSSWALK)[cohort=="LUAD" & target_gene=="STK11"]
if(nrow(cw)!=110L||anyDuplicated(cw$sample_id)) stop("STK11 crosswalk identity failure")
cw[,benchmark_group:=fifelse(canonical_lof==TRUE,"CANONICAL_LOF",fifelse(clean_wt==TRUE,"CLEAN_WT","EXCLUDED"))]
lof<-cw[benchmark_group=="CANONICAL_LOF"]
wt<-cw[benchmark_group=="CLEAN_WT"]
if(nrow(lof)!=10L||nrow(wt)!=88L) stop("frozen feasibility group count drift")
write_tsv(lof, file.path(OUT3,"CPTAC_CANONICAL_LOF_GROUPS.tsv"))
write_tsv(wt, file.path(OUT3,"CPTAC_CLEAN_WT_GROUPS.tsv"))

# RNA: exactly the frozen Core C network versions and run_ulm scoring.
expr<-read_symbol_matrix(RNA)
do<-fread(DOROTHEA); setnames(do,c("tf","confidence","target","mor"),c("source","confidence","target","mor")); do[,mor:=as.numeric(mor)];do<-do[confidence%chin%c("A","B")&is.finite(mor)]
pr<-fread(PROGENY);setnames(pr,c("gene","pathway","weight","p.value"),c("target","source","mor","p_value"));pr[,mor:=as.numeric(mor)];pr<-pr[is.finite(mor)]
tf_score<-as.data.table(run_ulm(expr,do[,.(source,target,mor)],.source="source",.target="target",.mor="mor",minsize=5L,na.rm=TRUE))
pr_score<-as.data.table(run_ulm(expr,pr[,.(source,target,mor)],.source="source",.target="target",.mor="mor",minsize=5L,na.rm=TRUE))
setnames(tf_score,c("source","condition","score"),c("feature","sample_id","score"),skip_absent=TRUE)
setnames(pr_score,c("source","condition","score"),c("feature","sample_id","score"),skip_absent=TRUE)
tf_score[,endpoint:=paste0("DOROTHEA_",feature)];pr_score[,endpoint:=paste0("PROGENY_",feature)]
selected<-rbind(tf_score[endpoint%chin%c("DOROTHEA_HIF1A","DOROTHEA_FOXO3"),.(sample_id,endpoint,score)],pr_score[endpoint%chin%c("PROGENY_PI3K","PROGENY_Hypoxia"),.(sample_id,endpoint,score)])
if(uniqueN(selected$endpoint)!=4L||anyDuplicated(selected,by=c("sample_id","endpoint"))) stop("frozen RNA axis feature failure")
selected[,z_score:=zscore(score),by=endpoint]
selected[,signed_z:=z_score*fifelse(endpoint=="DOROTHEA_FOXO3",-1,1)]
composite<-selected[,.(score=mean(signed_z),z_score=mean(signed_z),signed_z=mean(signed_z)),by=sample_id]
composite[,endpoint:="STK11_LKB1_AMPK_MTOR_METABOLIC"]
activity<-rbind(selected[,.(sample_id,endpoint,score,z_score,signed_z)],composite[,.(sample_id,endpoint,score,z_score,signed_z)])
activity<-merge(activity,cw[,.(sample_id,benchmark_group)],by="sample_id",all.x=TRUE)
fwrite(activity,file.path(OUT4,"CPTAC_RNA_ACTIVITY_MATRIX.tsv.gz"),sep="\t",quote=FALSE,na="NA",compress="gzip")

rna_defs<-data.table(endpoint=c("STK11_LKB1_AMPK_MTOR_METABOLIC","PROGENY_PI3K","PROGENY_Hypoxia","DOROTHEA_HIF1A","DOROTHEA_FOXO3"),expected_sign=c(1,1,1,1,-1),role=c("PRIMARY","COMPONENT","COMPONENT","COMPONENT","COMPONENT"))
rna_results<-rbindlist(lapply(seq_len(nrow(rna_defs)),function(i){d<-rna_defs[i];v<-activity[endpoint==d$endpoint];summarize_endpoint(v$z_score,v$benchmark_group,d$endpoint,"RNA",d$expected_sign,230828L+i,d$role)}))
write_tsv(rna_results,file.path(OUT4,"CPTAC_RNA_AXIS_RESULTS.tsv"))

# Proteomics: target protein plus the prespecified six-protein programme.
pmat<-read_symbol_matrix(PROTEIN)
samples<-colnames(pmat); pdt<-data.table(sample_id=samples)
target<-if("STK11"%in%rownames(pmat)) as.numeric(pmat["STK11",]) else rep(NA_real_,length(samples))
pdt[,STK11_TARGET_PROTEIN:=target]
panel<-c("FASN","ACACA","ACLY","SCD","RPS6KB1","EIF4EBP1")
present<-intersect(panel,rownames(pmat))
if(length(present)<4L) stop("fewer than four frozen proteomic programme components quantified")
panel_z<-t(apply(pmat[present,,drop=FALSE],1,zscore)); colnames(panel_z)<-samples
pdt[,STK11_MTORC1_LIPOGENIC_PROGRAMME:=apply(panel_z,2,function(v) if(sum(is.finite(v))>=4L) mean(v,na.rm=TRUE) else NA_real_)]
pdt<-merge(pdt,cw[,.(sample_id,benchmark_group)],by="sample_id",all.x=TRUE)
target_result<-summarize_endpoint(pdt$STK11_TARGET_PROTEIN,pdt$benchmark_group,"STK11_TARGET_PROTEIN","PROTEIN",-1,230840L,"PRIMARY")
program_result<-summarize_endpoint(pdt$STK11_MTORC1_LIPOGENIC_PROGRAMME,pdt$benchmark_group,"STK11_MTORC1_LIPOGENIC_PROGRAMME","PROTEIN",1,230841L,"PRIMARY")
target_result[,target_protein_status:=ifelse(n_lof>=8L&&n_wt>=30L,"QUANTIFIED","NOT_QUANTIFIED_OR_INSUFFICIENT")]
program_result[,components_frozen:=paste(panel,collapse=";")];program_result[,components_quantified:=paste(present,collapse=";")]
write_tsv(target_result,file.path(OUT5,"CPTAC_TARGET_PROTEIN_RESULTS.tsv"))
write_tsv(program_result,file.path(OUT5,"CPTAC_PROTEOMIC_PROGRAMME_RESULTS.tsv"))

rna_class<-rna_results[endpoint=="STK11_LKB1_AMPK_MTOR_METABOLIC",replication_class]
pclasses<-c(target_result$replication_class,program_result$replication_class)
integrated<-if(rna_class=="DIRECTIONALLY_REPLICATED"&&any(pclasses=="DIRECTIONALLY_REPLICATED")) "RNA_AND_PROTEIN_REPLICATED" else if(rna_class=="DIRECTIONALLY_REPLICATED") "RNA_REPLICATED_PROTEIN_UNRESOLVED" else if(rna_class=="OPPOSITE_DIRECTION"||all(pclasses=="OPPOSITE_DIRECTION")) "OPPOSITE_DIRECTION" else "NO_CLEAR_REPLICATION"
final_state<-if(integrated=="RNA_AND_PROTEIN_REPLICATED") "CPTAC_EXTERNAL_VALIDATION_SUCCESS" else if(rna_class%chin%c("DIRECTIONALLY_REPLICATED","PARTIALLY_REPLICATED")||any(pclasses%chin%c("DIRECTIONALLY_REPLICATED","PARTIALLY_REPLICATED"))) "CPTAC_EXTERNAL_VALIDATION_PARTIAL" else "FINAL_EXTENSION_NO_ADDITIONAL_SUPPORT"
integrated_tab<-data.table(cohort="LUAD",gene="STK11",canonical_lof_n=nrow(lof),clean_wt_n=nrow(wt),feasibility="GO_FULL_MULTIOMIC",rna_replication=rna_class,target_protein_replication=target_result$replication_class,proteomic_programme_replication=program_result$replication_class,integrated_class=integrated,final_extension_state=final_state,claim_ceiling="independent observational replication of canonical STK11-loss functional state; not validation of frozen9 or causality")
write_tsv(integrated_tab,file.path(OUT6,"CPTAC_INTEGRATED_GENE_RESULTS.tsv"))

writeLines(c("# CPTAC RNA replication report","",paste0("STK11 primary RNA axis: `",rna_class,"`."),"","The axis was frozen before outcome inspection and uses the same DoRothEA A/B, PROGENy top-100 and decoupleR run_ulm framework as Core C. Wilcoxon P is secondary; replication classification uses direction, bootstrap uncertainty and Cliff's delta."),file.path(OUT4,"CPTAC_RNA_REPLICATION_REPORT.md"))
writeLines(c("# CPTAC proteomic validation report","",paste0("STK11 target protein: `",target_result$replication_class,"`."),paste0("Prespecified mTORC1/lipogenic programme: `",program_result$replication_class,"`."),"","Total protein abundance is orthogonal observational support and is not phospho-activity or causal evidence. Phosphoproteomics was skipped as not required."),file.path(OUT5,"CPTAC_PROTEOMIC_VALIDATION_REPORT.md"))
writeLines(c("# CPTAC external validation report","",paste0("Integrated STK11 result: `",integrated,"`."),paste0("Final extension state: `",final_state,"`."),"","This independently assesses canonical STK11 loss biology in CPTAC-LUAD; it does not externally validate the nine TCGA R2 events."),file.path(OUT6,"CPTAC_EXTERNAL_VALIDATION_REPORT.md"))

audit<-data.table(metric=c("status","decoupleR_version","data.table_version","bootstrap_iterations","bootstrap_seed_base","rna_matrix_samples","protein_matrix_samples","RNA_axis_components","protein_programme_components_frozen","protein_programme_components_quantified","phospho_status","fallback_executed","outcome_scope"),value=c("PASS",as.character(packageVersion("decoupleR")),as.character(packageVersion("data.table")),B,230828,ncol(expr),ncol(pmat),4,6,length(present),"SKIPPED_NOT_REQUIRED","FALSE","STK11_LUAD_ONLY_AFTER_FEASIBILITY_GO"))
write_tsv(audit,file.path(OUT6,"CPTAC_VALIDATION_AUDIT.tsv"))
cat("CPTAC_STK11_FORMAL_STATUS=PASS\n")
cat("STK11_RNA_REPLICATION=",rna_class,"\n",sep="")
cat("STK11_TARGET_PROTEIN_REPLICATION=",target_result$replication_class,"\n",sep="")
cat("STK11_PROTEOMIC_PROGRAMME_REPLICATION=",program_result$replication_class,"\n",sep="")
cat("CPTAC_INTEGRATED_CLASS=",integrated,"\n",sep="")
cat("FINAL_EXTENSION_STATE=",final_state,"\n",sep="")
