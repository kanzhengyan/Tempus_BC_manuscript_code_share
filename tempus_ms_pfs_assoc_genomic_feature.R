#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 17, 2023
#
# script to calculate PFS associations from genomic features
# pre-treatment samples only
#
# tempus_ms_pfs_assoc_genomic_feature.R
#
# - parse a list of binary variables (e.g. genomic features)
# - calculate PFS associations for all features
# - create km plots for binary features
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2)
library(survival);
library(survminer);
source('Code/TEMPUS_BC.R');

infile <- paste(outdir, "tempus_ms_genomic_features.txt", sep='/')
outfile <- paste(outdir, "tempus_ms_pfs_assoc_genomic_features.txt", sep='/')
plotfile <- paste(outdir, 'tempus_ms_pfs_assoc_genomic_feature_kmplots.pdf', sep='/');

#--------------------------------
# (1) get tempus data (pfs)
#
#rdata = tempus_bc();
sdata = rdata$sdata;

# revise data
sdata$pfs_censored <- sdata$pfs_event;

#--------------------------------
# (2) parse genomic features
#
tmp <- read.table(infile, header=T, sep='\t');
mdata <- merge(sdata[,c('sample_id', 'pfs', 'pfs_censored', 'prepost')], 
               data.frame(sample_id=rownames(tmp), tmp), by='sample_id')
features <- colnames(tmp)

#----------------------------------------------
# (3) calculate PFS association statistics
#
sidx_pre <- (mdata$prepost == 'Pre');
foo <- mdata[sidx_pre,c('pfs', 'pfs_censored')];
data <- mdata[sidx_pre, features];

# calculate PFS association for each feature
res_binary <- apply(data, 2, function(x){
  
  # p-value log-rank test
  sdf <- survdiff(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1);
  
  # median survival
  sft <- survfit(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  n_wt <- summary(sft)$table[,"records"][1];
  n_alt <- summary(sft)$table[,"records"][2];
  mpfs_wt <- summary(sft)$table[,"median"][1];
  mpfs_alt <- summary(sft)$table[,"median"][2];
  
  # hazard ratio (HR) from CPH
  cph <- coxph(Surv(pfs, pfs_censored)~(x == 1), data=foo);
  hr <- round(coef(summary(cph))[,2],3);
  
  return(data.frame(pval=p.val, statistic=sdf$chisq, n_wt, n_alt, mpfs_wt, mpfs_alt, hr));}
)
res <- do.call("rbind", res_binary)
res$logp <- ifelse(res$hr > 1, -log10(res$pval), log10(res$pval))
res$type <- ifelse(grepl('MUT.', rownames(res)), 'MUT',
                   ifelse(grepl('AMP.', rownames(res)), 'AMP',
                          ifelse(grepl('DEL.', rownames(res)), 'DEL', '')));
res$gene <- gsub('DEL.', '', gsub('AMP.', '', gsub('MUT.', '', rownames(res))))

# FDR adjust for each alteration type
tmp1 <- res[res$type == 'MUT',];
tmp1$qval <- p.adjust(tmp1$pval, method='fdr')
tmp2 <- res[res$type == 'AMP',];
tmp2$qval <- p.adjust(tmp2$pval, method='fdr')
tmp3 <- res[res$type == 'DEL',];
tmp3$qval <- p.adjust(tmp3$pval, method='fdr')
tmp <- rbind(tmp1, tmp2, tmp3)
out <- tmp[,c('type', 'gene', 'pval', 'statistic', 'n_wt', 'n_alt', 
              'mpfs_wt', 'mpfs_alt', 'hr', 'logp', 'qval')];

#------------------------------
# write to output file
#
write.table(out, file=outfile, quote=F, sep="\t", row.names=F, col.names=T)

#----------------------------------------------
# (4) create KM plots for discrete features
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=8, pointsize=12);

# genomic features w/ significant PFS assoc.
features <- c('MUT.TP53', 'MUT.AKT1', 'MUT.NF1', 'MUT.BRCA2', 'MUT.PIK3CA', 'MUT.IRS2',
              'MUT.LRP1B', 'MUT.BCORL1', 'MUT.KDM6A', 'MUT.TAF1', 'MUT.PBRM1', 'MUT.CSF3R',
              'AMP.RAD51C', 'AMP.NTHL1', 'AMP.FGFR1', 'AMP.CCND1', 'AMP.MDM2', 'DEL.CDKN2A');
mdata_pre = mdata[sidx_pre,];

for (i in 1:length(features)) {
  feature = features[i];
  print(paste(i, feature, sep=': '));
  tmp <- mdata_pre[,feature];

  if (grepl('^MUT', feature)){
    vals <- ifelse(tmp == 1, 'MUT', 'WT');
    } else if(grepl('^AMP', feature)) {
      vals <- ifelse(tmp == 1, 'AMP', 'WT');
    } else if(grepl('^DEL', feature)){
      vals <- ifelse(tmp == 1, 'DEL', 'WT');
    }

  vardata = data.frame(sample_id=mdata_pre$sample_id, value=vals);
  foo <- merge(mdata_pre, vardata, by='sample_id', sort=F);
  fit <- survfit(as.formula(paste('Surv(pfs, pfs_censored)~', 'value')), data=foo);
  p <- ggsurvplot(fit, conf.int = FALSE, ggtheme = theme_classic(base_size = 19), palette=c('red', 'navyblue'), 
                  title=paste(feature, paste('(n=', sum(sidx_pre), ')', sep='')),
                  legend.labs=names(table(foo$value)), legend.title = '',
                  ylab = "PFS probability", xlab = "Time (days)", 
                  pval=F, risk.table=T, fontsize=7, font.tickslab=22, size=2, censor.size=10);
  
  # annotate p-value, HR from CPH analysis
  r <- res[rownames(res) == feature,]
  label1 <- paste('p=', format(r$pval, digits=4), sep='');
  label2 <- paste('HR=', format(r$hr, digits=4), sep='');
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 1, label = label1, size = 8) # x and y coordinates of the text
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 0.9, label = label2, size = 8) # x and y coordinates of the text
  
  print(p);
}

dev.off();