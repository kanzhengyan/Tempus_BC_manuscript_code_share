#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# Nov. 25, 2023
# Sep. 2, 2024
#
# script to characterize gene expression and gene signatures
# using multi-panel boxplots
#
setwd('../Fig3.FigS9.cse_gsva_DE_analysis')
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(RColorBrewer);

plotfile <- paste(outdir, 'Fig-3.cse_gsva_waterfall_analysis.pdf', sep='/');

#--------------------------------
# Tempus BC R data object
#
#rdata = tempus_bc();
#sdata = rdata$sdata;

#---------------------------------------------
# (1) get sample metadata
#
mfdatafile <- paste(datadir, "Tempus_molecular_features_extended.txt", sep='/');
mdata <- read.table(mfdatafile, sep="\t", header=T);

mdata$prepost = factor(mdata$pre_post, levels=c('Pre', 'Post'))
mdata$icluster2 = ifelse(mdata$icluster == 'IC1', 'IC1', 'IC2-5');
mdata <- subset(mdata, select=c('sample_id', 'patient_id', 'subcohorts', 'prepost',
                                'subtype_pam50', 'icluster', 'icluster2'));

#------------------------------------------------------
# (2) add CSE expression & gsva data
#

# load CSE GSVA data
gsvafile <- paste(datadir, 'tempus_ms_supp_data.hallmark_cse_gsva.txt', sep='/');
gsvadata <- read.table(gsvafile, sep='\t', header=T)
rownames(gsvadata) <- gsvadata[,1]
gsva <- t(gsvadata[,-1]);

# add gene signature GSVA data
rownames(gsva) <- gsub('EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT', rownames(gsva))
df2 <- data.frame(sample_id=colnames(gsva), t(gsva));
mdata <- merge(mdata, df2, by='sample_id', all.x=T, sort=F);
gsvanames <- rownames(gsva)

# create waterfall plot - signatures ranked by log10p
pdf(file=plotfile, onefile=T, height=8, width=11, pointsize=12);

#----------------------------------------
# (3) Compare mean GSVA score
# IC1 vs. IC2-5, for each signature
#
df <- NULL;
for (i in 1:length(gsvanames)){
  yvar <- gsvanames[i];
  a <- mdata[mdata$icluster2 == 'IC1', yvar];
  b <- mdata[mdata$icluster2 == 'IC2-5', yvar];
  res <- wilcox.test(a, b);
  df <- rbind(df, data.frame(yvar, variable="icluster2", 
                             a=mean(a, na.rm=T), b=mean(b, na.rm=T), pval=res$p.value));
}
df$logp <- ifelse(df$a > df$b, -log10(df$pval), log10(df$pval));
df$qval <- p.adjust(df$pval, method='BH')
rownames(df) <- gsvanames;

# order by difference in mean scores (log10P)
df <- df[order(df$logp, decreasing=T),]

# color by significance & direction
cols <- brewer.pal(4, 'Set1');
collabs <- ifelse (df$logp > 0 & df$pval < 0.05, cols[1], 
                   ifelse (df$logp < 0 & df$pval < 0.05, cols[2], 'grey'))

par(mar=c(20,8,4,4))
title <- 'IC1 vs. IC2-5';
barplot(df$logp, beside=T, xlab='', ylab='Signed Log10(p)\nIC1 vs.IC2-5', las=2, col=collabs,
        names.arg=rownames(df), cex.axis=1.5, cex.names=1, cex.lab=1.5, main=title);
abline(h=log10(0.05), lty=2);
abline(h=-log10(0.05), lty=2);

#----------------------------------------
# (4) Compare mean GSVA score
# Basal vs. other subtypes, for each signature
#
df <- NULL;
for (i in 1:length(gsvanames)){
  yvar <- gsvanames[i];
  a <- mdata[mdata$subtype_pam50 == 'Basal', yvar];
  b <- mdata[mdata$subtype_pam50 %in% c('Her2', 'LumA', 'LumB'), yvar];
  res <- wilcox.test(a, b);
  df <- rbind(df, data.frame(yvar, variable="pam50", 
                             a=mean(a, na.rm=T), b=mean(b, na.rm=T), pval=res$p.value));
}
df$logp <- ifelse(df$a > df$b, -log10(df$pval), log10(df$pval));
df$qval <- p.adjust(df$pval, method='BH')
rownames(df) <- gsvanames;

# order by difference in mean scores (log10P)
df <- df[order(df$logp, decreasing=T),]

# color by significance & direction
cols <- brewer.pal(4, 'Set1');
collabs <- ifelse (df$logp > 0 & df$pval < 0.05, cols[1], 
                   ifelse (df$logp < 0 & df$pval < 0.05, cols[2], 'grey'))

par(mar=c(20,8,4,4))
title <- 'Basal vs. other PAM50 subtypes';
barplot(df$logp, beside=T, xlab='', ylab='Signed Log10(p)\nBasal vs. Others', las=2, col=collabs,
        names.arg=rownames(df), cex.axis=1.5, cex.names=1, cex.lab=1.5, main=title);
abline(h=log10(0.05), lty=2);
abline(h=-log10(0.05), lty=2);

dev.off()

