#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 23, 2023
#
# script to calculate PFS associations from gene signatures
# pre-treatment samples only
#
# tempus_ms_pfs_assoc_gsva.R
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2)
library(survival);
library(survminer);
source('Code/TEMPUS_BC.R');

outfile <- paste(outdir, "tempus_ms_pfs_assoc_gsva.txt", sep='/')
plotfile <- paste(outdir, 'tempus_ms_pfs_assoc_gsva_kmplots.pdf', sep='/');

#--------------------------------
# (1) get tempus data (pfs)
#
#rdata = tempus_bc();
sdata = rdata$sdata;

# revise data
sdata$pfs_censored <- sdata$pfs_event;

#------------------------------------
# (2) load gene signatures (gsva)
#
gsvaobj <- tempus_bc_gsva(rdata);
gsvadata <- gsvaobj$scores;
gidx <- grepl('h.HALLMARK_', rownames(gsvadata))

# align with sample annot
midx <- match(sdata$sample_id, colnames(gsvadata));
gsva <- t(gsvadata[gidx,midx]);
colnames(gsva) <- gsub('h.HALLMARK_', '', colnames(gsva));

mdata <- merge(sdata[,c('sample_id', 'pfs', 'pfs_censored', 'prepost')], 
               data.frame(sample_id=rownames(gsva), gsva), by='sample_id')
features <- colnames(gsva)

#----------------------------------------------
# (3) calculate PFS association statistics
#
sidx_pre <- (mdata$prepost == 'Pre');
foo <- mdata[sidx_pre,c('pfs', 'pfs_censored')];
data <- mdata[sidx_pre, features];

# convert continuous variable to binary (median-split)
tmp <- apply(data, 2, function(x) return(x > quantile(x, probs=c(0.5), na.rm=T)))
medians <- apply(data, 2, function(x) return(round(quantile(x, probs=c(0.5), na.rm=T), 4)));
data2 <- ifelse(tmp, 1, 0)

# calculate PFS association for each feature
res_binary <- apply(data2, 2, function(x){
  
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
df <- data.frame(variable=rownames(res), res);

#------------------------------
# write to output file
#
write.table(df, file=outfile, quote=F, sep="\t", row.names=F, col.names=T)

#----------------------------------------------
# (4) create KM plots for discrete features
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=8, pointsize=12);

# genomic features w/ significant PFS assoc.
mdata_pre = mdata[sidx_pre,];
pcut <- 0.05;
features_sig <- df[order(df$pval) & df$pval < pcut,]$variable;

# add custom signatures
features_sig <- c('E2F_TARGETS', 'MYC_TARGETS_V1', 'G2M_CHECKPOINT', 'PI3K_AKT_MTOR_SIGNALING',
                 features_sig)

for (i in 1:length(features_sig)) {
  feature = features_sig[i];
  print(paste(i, feature, sep=': '));
  tmp <- data2[,feature];
  vals <- ifelse(tmp == 1, paste('High (>', medians[feature], ')', sep=''), 
                 paste('Low (<=', medians[feature], ')', sep=''));

  vardata = data.frame(sample_id=mdata_pre$sample_id, value=vals);
  foo <- merge(mdata_pre, vardata, by='sample_id', sort=F);
  fit <- survfit(as.formula(paste('Surv(pfs, pfs_censored)~', 'value')), data=foo);
  p <- ggsurvplot(fit, conf.int = FALSE, ggtheme = theme_classic(base_size = 19), palette=c('red', 'navyblue'), 
                  title=feature, legend.labs=names(table(foo$value)), legend.title = '',
                  ylab = "PFS probability", xlab = "Time (days)", 
                  pval=F, risk.table=T, fontsize=7, font.tickslab=22, size=2, censor.size=10);
  
  # annotate p-value, HR from CPH analysis
  r <- df[df$variable == feature,]
  label1 <- paste('p=', format(r$pval, digits=4), sep='');
  label2 <- paste('HR=', format(r$hr, digits=4), sep='');
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 1, label = label1, size = 8) # x and y coordinates of the text
  p$plot  <- p$plot + ggplot2::annotate("text", x = 1500, y = 0.9, label = label2, size = 8) # x and y coordinates of the text
  
  print(p);
}

dev.off();
