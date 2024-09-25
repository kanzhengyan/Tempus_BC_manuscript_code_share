#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# August 15, 2023
# Sept. 21, 2023
#
# script to characterize predicted perturbation response scores
#
setwd('../Fig5.ENRS_characterize');
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library('reshape2');
library('RColorBrewer');

#--------------------------------
# Tempus BC R data object
#
#rdata = tempus_bc();
#sdata = rdata$sdata;

plotfile <- paste(outdir, 'Fig-5.ENRS_score_characterize.pdf', sep='/');

#------------------------------------------------------
# (1) get predicted perturbation response scores
# distribution & correlation patterns
#
mfdatafile <- paste(datadir, "Tempus_molecular_features_extended.txt", sep='/');
mdata <- read.table(mfdatafile, sep="\t", header=T);
mdata$prepost <- mdata$pre_post;
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;

# load response scores for analysis
enrsfile <- paste(datadir, 'tempus_ms_supp_data.enrs_score.txt', sep='/');
enrsdata <- read.table(enrsfile, sep='\t', header=T)
scores = colnames(enrsdata)[grepl('^ENRS', colnames(enrsdata))]
mdata <- merge(mdata, enrsdata, by='sample_id', all.x=T, sort=F)

#============================================
# ENRS score characterization analysis
#
pdf(file=plotfile, onefile=T, height=8, width=12, pointsize=12);

# color labels
annot_colors = list(
  icluster =  c(brewer.pal(4, 'Set1'), 'grey'),
  prepost = c('Post'='grey20', 'Pre'='lightgrey'),
  pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  histology = rainbow(3),
  organ = rainbow(4)
)

# reorder levels
mdata$prepost = factor(mdata$prepost, levels=c('Pre', 'Post'))  

#---------------------------------------------------------
# (1) compare variable vs. key categories 
# Pre/Post, Pre/During/Post, Tissue, PAM50, IC, 
# MUT, AMP
#
plist <- plist2 <- list();
c = 0;
categories <- c('prepost', 'subtype_pam50', 'icluster');

for (i in 1:length(categories)){
  category <- categories[i];

  # append data for multiple variables
  for (j in 1:length(scores)){
    score <- scores[j];
    df <- data.frame(mdata, score=mdata[,score], varlabel=gsub('ENRS_CERES_', '', score));
    if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
  }
  
  # reorder levels
  foo$varlabel = factor(foo$varlabel, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))  
  
  if (category == 'prepost') {
    comps <- list(c("Pre", "Post"));
    orders <- c('Pre', 'Post');
    colors <- annot_colors$prepost;
  } else if (category == 'subtype_pam50') {
    comps <- list(c('LumA', 'LumB'), c('LumA', 'Her2'), c('LumA', 'Basal'));
    orders <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
    colors <- annot_colors$pam50;
    foo <- foo[foo$subtype != 'Normal',];
  } else if (category == 'icluster') {
    comps <- list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'), c('IC1', 'IC5'));
    orders <- c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
    colors <- annot_colors$icluster;
    foo <- foo[foo$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4', 'IC5'),]
  }

  plist[[i]] = ggboxplot(foo, category, 'score', xlab='', ylab='Dependency Score', outlier.shape=NA,
                         facet.by='varlabel', scales='free_y', ncol=length(scores),
                         fill=category, palette=colors,  order=orders) +
    stat_compare_means(comparison=comps, size=5, label='p.signif') +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    theme(text = element_text(size = 12, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2)
}
print(ggarrange(plotlist=plist, nrow=2, ncol=1));
 
#---------------------------------------------------------------------------
# (2) Longitudinal variable changes Pre vs. Post (paired BL/PT)
#

# data for only sample pairs
mdata_paired = mdata[mdata$subcohorts == '1-Paired',]

# append data for multiple variables
for (j in 1:length(scores)){
  score <- scores[j];
  df <- data.frame(mdata_paired, score=mdata_paired[,score], varlabel=gsub('ENRS_CERES_', '', score));
  if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
}
foo$varlabel = factor(foo$varlabel, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))  

p <- ggpaired(foo, x='prepost', y='score', title='', fill='prepost', id='patient_id', xlab='', ylab='Dependency Score',
               facet.by='varlabel', scales='free_y', ncol=length(scores),
               line.color='gray', line.size=.4, palette=annot_colors$prepost, order=c('Pre', 'Post')) +
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  stat_compare_means(paired=TRUE, size=5, label='p.format') +
  theme(text = element_text(size = 12, face='bold'), legend.position = "none");

print(ggarrange(p, nrow=2, ncol=1));

#------------------------------------------------------------------
# (3) compare variable vs. key categories 
# bar plot of statistical significance of
# comparing mean dependency score for MUT vs. WT of each gene
#
genomic_muts = c('MUT.ESR1', 'MUT.RB1', 'MUT.TP53');
mm <- melt(mdata[,c('sample_id', genomic_muts)], id.vars=1);
mm$value <- ifelse(mm$value == 1, 'MUT', 'WT');
mm$variable <- gsub('MUT.', '', mm$variable);

df <- NULL;
for (i in 1:length(scores)){
  yvar <- scores[i];
  foo <- merge(mm, mdata[,c('sample_id', yvar)], by='sample_id')

  genes <- levels(as.factor(foo$variable));
  for (j in 1:length(genes)){
    gene <- genes[j];
    a <- foo[foo$variable == gene & foo$value == 'MUT', yvar];
    b <- foo[foo$variable == gene & foo$value == 'WT', yvar];
    res <- wilcox.test(a, b);
    df <- rbind(df, data.frame(yvar, gene, a=mean(a), b=mean(b), pval=res$p.value));
  }
}
df$logp <- ifelse(df$a > df$b, -log10(df$pval), log10(df$pval));
df$ratio <- ifelse(df$a > df$b, abs(log2(df$a/df$b)), -abs(log2(df$a/df$b)))
df$delta <- df$a - df$b;

par(mar=c(12,12,4,4));
df$yvar <- gsub('ENRS_CERES_', '', df$yvar);
df$yvar <- factor(df$yvar, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))

# signed log(p)
barplot(logp~yvar+gene, beside=T, xlab='', ylab='Signed Log10(p)\nMUT vs.WT', data=df, las=1, col=brewer.pal(4, 'Set1'), 
        cex.axis=2, cex.names=2, cex.lab=1.5, ylim=c(-15,15));
abline(h=log10(0.05), lty=2);
abline(h=-log10(0.05), lty=2);
legend("topleft", legend=levels(as.factor(df$yvar)), fill=brewer.pal(4, 'Set1'), 
       title="Gene Dependency", ncol=2, cex=1.2, bg='white');

# log2(ratio)
barplot(ratio~yvar+gene, beside=T, xlab='', ylab='Dependency Ratio (Log2)\nMUT vs.WT', data=df, las=1, col=brewer.pal(4, 'Set1'), 
        cex.axis=1.5, cex.names=2, cex.lab=1.5, ylim=c(-.2, .2));
legend("topleft", legend=levels(as.factor(df$yvar)), fill=brewer.pal(4, 'Set1'), 
       title="Dependency Ratio", ncol=2, cex=1.2, bg='white');

#-----------------------------------------------
# (6) comparing mean dependency score 
# Pre vs. Post, for each gene
#
df <- NULL;
for (i in 1:length(scores)){
  yvar <- scores[i];
  a <- mdata[mdata$pre_post == 'Pre', yvar];
  b <- mdata[mdata$pre_post == 'Post', yvar];
  res <- wilcox.test(a, b);
  df <- rbind(df, data.frame(yvar, variable="PrePost", a=mean(a), b=mean(b), pval=res$p.value));
}
df$logp <- ifelse(df$a < df$b, -log10(df$pval), log10(df$pval));

par(mar=c(2,7,2,5), mfrow=c(2,3));
df$yvar <- gsub('ENRS_CERES_', '', df$yvar);
df$yvar <- factor(df$yvar, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))
barplot(logp~yvar, beside=T, xlab='', ylab='Signed Log10(p)\nPre vs.Post', data=df, las=1, col=brewer.pal(4, 'Set1'), 
        xaxt='n', cex.axis=1.5, cex.names=1.4, cex.lab=1.5, ylim=c(-15, 15));
abline(h=log10(0.05), lty=2);
abline(h=-log10(0.05), lty=2);
legend("topright", legend=levels(as.factor(df$yvar)), fill=brewer.pal(4, 'Set1'), 
       title="Gene Dependency", ncol=2, cex=1, bg='white');

#-----------------------------------------------
# (7) comparing mean dependency score 
# LumB: IC1 vs. IC2-4, for each gene
#
df <- NULL;
for (i in 1:length(scores)){
  yvar <- scores[i];
  a <- mdata[mdata$icluster %in% c('IC2', 'IC3', 'IC4', 'IC5') & mdata$subtype_pam50 %in% c('LumA', 'LumB'), yvar];
  b <- mdata[mdata$icluster == 'IC1' & mdata$subtype_pam50 %in% c('LumA', 'LumB'), yvar];
  res <- wilcox.test(a, b);
  df <- rbind(df, data.frame(yvar, variable="PAM50-IC", a=mean(a, na.rm=T), b=mean(b, na.rm=T), pval=res$p.value));
}
df$logp <- ifelse(df$a < df$b, -log10(df$pval), log10(df$pval));

df$yvar <- gsub('ENRS_CERES_', '', df$yvar);
df$yvar <- factor(df$yvar, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))
barplot(logp~yvar, beside=T, xlab='', ylab='Signed Log10(p)\nLum:IC1 vs. IC2-5', data=df, las=1, col=brewer.pal(4, 'Set1'),
        xaxt="n", cex.axis=1.5, cex.names=1.4, cex.lab=1.5, ylim=c(-8, 8));
abline(h=log10(0.05), lty=2);
abline(h=-log10(0.05), lty=2);
legend("topright", legend=levels(as.factor(df$yvar)), fill=brewer.pal(4, 'Set1'),
       title="Gene Dependency", ncol=2, cex=1, bg='white');

dev.off();

