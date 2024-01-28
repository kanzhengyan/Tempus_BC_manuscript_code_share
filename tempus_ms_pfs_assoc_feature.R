#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 10, 2023
# May 29, 2023
#
# script to perform PFS association analysis &
# create KM plots from molecular features
# pre-treatment samples only
#
# tempus_ms_pfs_assoc_feature.R
#
# - extract a list of molecular features
# - convert continuous to binary features (median-split)
# - calculate PFS associations for all features
# - create km plots for binary features
# - create km plots for (multi-) categorical features
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2)
library(survival);
library(survminer);
source('Code/TEMPUS_BC.R');

outfile <- paste(outdir, "tempus_pfs_assoc_results_features.txt", sep='/')
plotfile <- paste(outdir, 'tempus_pfs_assoc_feature_kmplots.pdf', sep='/');

#----------------------------------------
# (1) load Rdata and molecular data
#
rdata = tempus_bc();
sdata = rdata$sdata;

# sample annotation with covariates
annotation <- sdata[, c('sample_id', 'patient_id', 'prepost', 'tumor_purity', 
                        'organ', 'proliferative_index', 'subcohorts')];
annotation$time <- factor(annotation$prepost, levels=c("Pre", "Post"))

# get molecular features and feature metadata
molfeatfile <- paste(datadir, 'Tempus_molecular_features.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');
metadatafile <- paste(datadir, 'Tempus_molecular_features_metadata.txt', sep='/');
meta <- read.table(metadatafile, header=T, sep='\t');

# align with sample annot
midx <- match(sdata$sample_id, mdata$sample_id)
mdata <- mdata[midx,];

# revise data
mdata$pfs <- mdata$pfs1;
mdata$pfs_censored <- mdata$pfs1_censored;
mdata <- mdata[,colnames(mdata) != 'msi_status']

#----------------------------------------------
# (2) calculate PFS association statistics
#

# extract binary & continuous variables
mode = 'all'
if (mode == 'all'){
  idx1 <- colnames(mdata) %in% meta$variable[meta$variable_type == 'binary' & 
                                               meta$feature_type %in% c('molecular', 'gene expression', 'genomic')]
  idx2 <- colnames(mdata) %in% meta$variable[meta$variable_type == 'numerical' & 
                                               meta$feature_type %in% c('molecular', 'gene expression', 'genomic')]
} else {
  idx1 <- colnames(mdata) %in% meta$variable[meta$variable_type == 'binary' & meta$feature_type %in% c('molecular')]
  idx2 <- colnames(mdata) %in% meta$variable[meta$variable_type == 'numerical' & meta$feature_type %in% c('molecular')]
}
features_binary <- colnames(mdata)[idx1]
features_continuous <- colnames(mdata)[idx2]

sidx_pre <- (mdata$pre_post == 'Pre');
foo <- mdata[sidx_pre, c('pfs', 'pfs_censored')];

# binary variables
data1 <- mdata[sidx_pre, features_binary];
data1$brca_pathogenic_mutation <- ifelse(data1$brca_pathogenic_mutation == 'Yes', 1, 0)
data1$bc_susceptibility_gene_pathogenic_mutation <- ifelse(data1$bc_susceptibility_gene_pathogenic_mutation == 'Yes', 1, 0)
# remove variables missing one outcome
data1 <- data1[, apply(data1, 2, sum, na.rm=T) >= 2];

# convert continuous variable to binary (median-split)
data2 <- mdata[sidx_pre, features_continuous];
tmp <- apply(data2, 2, function(x) return(x > quantile(x, probs=c(0.5), na.rm=T)))
medians <- apply(data2, 2, function(x) return(round(quantile(x, probs=c(0.5), na.rm=T), 4)));
data2 <- ifelse(tmp, 1, 0)
# remove variables missing one outcome
data2 <- data2[, apply(data2, 2, sum, na.rm=T) >= 2];

# combine binary & continous variables
data <- cbind(data1, data2)
rownames(data) <- mdata[sidx_pre,]$sample_id;
vartypes <- c(rep('binary', dim(data1)[2]),
              rep('continuous', dim(data2)[2]));
names(vartypes) <- colnames(data)

# calculate PFS association for each feature
res_binary <- apply(data, 2, function(x){
  
  # p-value log-rank test
  sdf <- survdiff(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1);
  
  # median survival
  sft <- survfit(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  n_dn <- summary(sft)$table[,"records"][1];
  n_up <- summary(sft)$table[,"records"][2];
  mpfs_dn <- summary(sft)$table[,"median"][1];
  mpfs_up <- summary(sft)$table[,"median"][2];
  
  # hazard ratio (HR) from CPH
  cph <- coxph(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  hr <- round(coef(summary(cph))[,2],3);
  
  return(data.frame(pval=p.val, statistic=sdf$chisq, n_dn, n_up, mpfs_dn, mpfs_up, hr));}
)
res <- do.call("rbind", res_binary)
res$qval <- p.adjust(res$pval)
res$logp <- ifelse(res$hr > 1, -log10(res$pval), log10(res$pval))
out <- data.frame(feature=rownames(res), res)

#------------------------------
# write to output file
#
write.table(out, file=outfile, quote=F, sep="\t", row.names=F, col.names=T)

#------------------------------------------------------
# (3) create KM plots for discrete features
#

#------------------------------------------------------
# Create KM plots for selected features
# binary: genomic alteration, PAM50 subtype
# continuous: gene expression, gene signature, ...
# PAM50 scores, NMF factors
#
if (mode == 'all') {

  # select from all features
  features_binary = c('MUT.TP53', 'MUT.RB1', 'MUT.PIK3CA', 'MUT.AKT1', 'MUT.ESR1',
                    'brca_pathogenic_mutation', 'MUT.PTEN', 'MUT.CDKN2A',
                    'AMP.CCND1', 'AMP.FGFR1', 'AMP.MDM2');
  features_continuous = c('proliferative_index', 'cyt_score', 'hrd_index', 'mutsig_13_sigma',
                        'pam50_cor_basal', 'pam50_cor_her2', 'pam50_cor_luma', 'pam50_cor_lumb', 'pam50_cor_normal',
                        'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT', 'Paloma3_F7_IFNG_response', 
                        'h.HALLMARK_MYC_TARGETS_V1', 'h.HALLMARK_ESTROGEN_RESPONSE_EARLY', 'h.HALLMARK_E2F_TARGETS', 
                        'c2_cgp.SMID_BREAST_CANCER_ERBB2_UP', 'h.HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'h.HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                        'EXPR.CCNE1', 'EXPR.ESR1', 'EXPR.PGR');
  features <- c(features_binary, features_continuous);
} else {

  # select only molecular features (analytical)
  features = c('brca_pathogenic_mutation', 'proliferative_index', 'cyt_score', 'mutsig_13_sigma',
                    'pam50_cor_basal', 'pam50_cor_her2', 'pam50_cor_luma', 'pam50_cor_lumb', 
                    'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 
                    'Paloma3_F8_EMT', 'Paloma3_F7_IFNG_response');
}
mdata_pre = mdata[sidx_pre, c('sample_id', 'pfs', 'pfs_censored')];

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=8, pointsize=12);

for (i in 1:length(features)) {
  feature = features[i];
  print(paste(i, feature, sep=': '));
  
  #----------------------------------------
  # (a) Pre survival analysis (KM plot)
  #
  tmp <- data[,feature];
  if (vartypes[feature] == 'continuous'){
    #vals <- ifelse(tmp == 1, 'High', 'Low');    
    vals <- ifelse(tmp == 1, paste('High (>', medians[feature], ')', sep=''), 
                   paste('Low (<=', medians[feature], ')', sep=''));
  } else{
    if (grepl('^MUT', feature) | grepl('mutation', feature)){
      vals <- ifelse(tmp == 1, 'MUT', 'WT');
    } else if(grepl('^AMP', feature)) {
      vals <- ifelse(tmp == 1, 'AMP', 'WT');
    } else if(grepl('^DEL', feature)){
      vals <- ifelse(tmp == 1, 'DEL', 'WT');
    }
  }
  vardata = data.frame(sample_id=rownames(data), value=vals);
  foo <- merge(mdata_pre, vardata, by='sample_id', sort=F);
  fit <- survfit(as.formula(paste('Surv(pfs, pfs_censored)~', 'value')), data=foo);
  title <- paste(feature, paste('(n=', sum(sidx_pre), ')', sep=''));
  title <- gsub('c2_cgp.', '', gsub('h.HALLMARK_', '', title));

  p <- ggsurvplot(fit, conf.int = FALSE, ggtheme = theme_classic(base_size = 18), palette=c('red', 'navyblue'), 
                   title=title,
                   legend.labs=names(table(foo$value)), legend.title = '',
                   ylab = "PFS probability", xlab = "Time (days)", 
                  pval=F, risk.table=T, fontsize=7, font.tickslab=22, size=2, censor.size=10);
  
  # annotate p-value, HR from CPH analysis
  label1 <- paste('p=', format(res[feature,]$pval, digits=4), sep='');
  label2 <- paste('HR=', format(res[feature,]$hr, digits=4), sep='');
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 1, label = label1, size = 8) # x and y coordinates of the text
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 0.9, label = label2, size = 8) # x and y coordinates of the text
    
  print(p);
}

#------------------------------------
# multi-categorical feature
#
feature <- 'subtype_pam50';
print(paste(i, feature, sep=': '));

# remove basal (n=1) and normal-like (n=1)
sidx_pam50 <- mdata$subtype_pam50 %in% c('LumA', 'LumB', 'Her2');
foo <- mdata[sidx_pre & sidx_pam50, c('sample_id', 'pfs', 'pfs_censored', feature)]
fit <- survfit(as.formula(paste('Surv(pfs, pfs_censored)~', feature)), data=foo);
p <- ggsurvplot(fit, conf.int = FALSE, ggtheme = theme_classic(base_size = 19), palette='lancet', 
                title=paste(feature, paste('(n=', sum(sidx_pre), ')', sep='')),
                legend.labs=names(table(foo[,feature])), legend.title = '',
                ylab = "PFS probability", xlab = "Time (days)", 
                pval=F, risk.table=T, fontsize=6, font.tickslab=22, size=2, censor.size=10);
p$plot <- p$plot + annotate("text", x = 1500, y = 1, label=surv_pvalue(fit)$pval.txt, size=8)
print(p)


dev.off();