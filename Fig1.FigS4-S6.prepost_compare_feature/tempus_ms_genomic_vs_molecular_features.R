#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 10, 2023
# Sept. 21, 2024
#
# script to create boxplots to show association
# for genomic alterations (discrete) vs. molecular features (continuous)
#
setwd('../Fig1.FigS4-S6.prepost_compare_feature')
datadir = 'Input';
outdir = 'Output.Manuscript';

library(reshape2)
library(ggpubr)

plotfile <- paste(outdir, 'Fig-S2.genomic_vs_molecular_features.pdf', sep='/');

#---------------------------------------------------------
# (1) get molecular features and feature metadata
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;
mdata$prepost <- mdata$pre_post;
mdata$tumor_purity <- mdata$tumor_purity_pathologist;
colnames(mdata)[grepl('h.HALLMARK', colnames(mdata))] <- gsub('h.HALLMARK_', '', colnames(mdata)[grepl('h.HALLMARK_', colnames(mdata))])
colnames(mdata)[grepl('c2_cgp.', colnames(mdata))] <- gsub('c2_cgp.', '', colnames(mdata)[grepl('c2_cgp.', colnames(mdata))])

#----------------------------------------------
# (2) create boxplots to show association
# mutations vs. expression features
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=10, pointsize=12);

genomic_muts = c('MUT.ESR1', 'MUT.RB1', 'MUT.TP53');
features = c('proliferative_index', 'ESTROGEN_RESPONSE_EARLY', 'EXPR.CCNE1');
colors_alter = c('skyblue2', 'firebrick3')

mm <- melt(mdata[,c('sample_id', genomic_muts)], id.vars=1)
mm$value <- ifelse(mm$value == 1, 'MUT', 'WT')

plist <- list();
comps <- list(c('WT', 'MUT'));
for (i in 1:length(features)){
  yvar <- features[i];
  foo <- merge(mm, mdata[,c('sample_id', yvar)], by='sample_id')
  plist[[i]] <- ggboxplot(foo, 'value', yvar, xlab='', ylab=yvar, title='', outlier.shape=NA, remove=NA,
                          fill='value', palette=colors_alter, facet.by='variable') + 
    stat_compare_means(size=4, label='p.format', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2);
}

ggarrange(plotlist=plist, nrow=2, ncol=2, legend=F);

dev.off();


