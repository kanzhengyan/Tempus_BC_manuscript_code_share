#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# August 15, 2023
#
# script to characterize predicted perturbation response scores
#
# tempus_ms_enrs_score_analysis.R > tempus_ms_enrs_score_analysis.pdf
#
# - extract a list of predicted response scores, plot
#   (1) score vs. Pre/Post time (boxplot, p-value)
#   (2) score vs. Pre/Post time (paired boxplot)
#   (3) score vs. Key categories (IC, IC-PAM50)
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library('reshape2');
library(RColorBrewer);
source('Code/TEMPUS_BC.R');

#--------------------------------
# Tempus BC R data object
#
rdata = tempus_bc();
sdata = rdata$sdata;

plotfile <- paste(outdir, 'tempus_ms_enrs_score_analysis.pdf', sep='/');

#------------------------------------------------------
# (1) get predicted perturbation response scores
# distribution & correlation patterns
#
mfdatafile <- paste(datadir, "Tempus_molecular_features_extended_v3.txt", sep='/');
mfdata <- read.table(mfdatafile, sep="\t", header=T);

# merge sample annotation w/ molecular features
mdata <- merge(mfdata, sdata[,c('sample_id', 'time', 'organ', 'subcohorts')], 
               by='sample_id', all.x=T, sort=F)
mdata$prepost <- mdata$pre_post;
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;

# select response scores for analysis  
scores = c('ENRS_CERES_ESR1', 'ENRS_CERES_CDK4', 'ENRS_CERES_CDK6', 'ENRS_CERES_CDK2');

pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);

#============================================
# function: score characterization analysis
#

# color labels
annot_colors = list(
  icluster =  brewer.pal(4, 'Set1'),
  prepost = c('Post'='grey20', 'Pre'='lightgrey'),
  time = c('BL'='white', 'OT'='lightgrey', 'PD'='grey3'),
  mut = c('0'='deepskyblue2', '1'='firebrick2'),
  amp = c('0'='deepskyblue2', '1'='orange3'),
  del = c('0'='deepskyblue2', '1'='darkblue'),
  pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  ic_pam50 = c('IC1:Basal'='Purple',  'IC1:Her2'='green', 'IC1:LumB'='orange', 
               'IC2:LumB'='orange', 'IC3:LumA'='yellow', 'IC3:LumB'='orange', 'IC4:LumA'='yellow'),
  brca_status = c('0'='deepskyblue2', '1'='firebrick2'),
  mutsig_status = rainbow(6),
  histology = rainbow(3),
  organ = rainbow(4)
)
sidx_pfs = grepl('PFS', mdata$pfs_status);
sidx_lum = mdata$subtype_pam50 %in% c('LumA', 'LumB');
  
# reorder levels
mdata$prepost = factor(mdata$prepost, levels=c('Pre', 'Post'))  

#-------------------------------------------
# add icluster+pam50 subgroups
#
groups <- paste(mdata$icluster, mdata$subtype_pam50, sep=':')
selected_groups <- names(table(groups))[table(groups) >= 7];
idx2 <- groups %in% selected_groups & mdata$icluster != 'IC5' & !is.na(mdata$icluster)
mdata$ic_pam50 <- groups;
mdata$ic_pam50[!idx2] <- NA;

#---------------------------------------------------------
# (1) compare variable vs. key categories 
# Pre/Post, Pre/During/Post, Tissue, PAM50, IC, 
# MUT, AMP
#
plist <- plist2 <- list();
c = 0;
categories <- c('prepost', 'subtype_pam50', 'icluster', 'time');

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
  } else if (category == 'time'){
    comps <- list(c("BL", "OT"), c("OT", "PD"), c("BL", "PD"));
    orders <- c('BL', 'OT', 'PD');
    colors <- annot_colors$time;
  } else if (category == 'subtype_pam50') {
    comps <- list(c('LumA', 'LumB'), c('LumA', 'Her2'), c('LumA', 'Basal'));
    orders <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
    colors <- annot_colors$pam50;
    foo <- foo[foo$subtype != 'Normal',];
  } else if (category == 'icluster') {
    comps <- list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'));
    orders <- c('IC1', 'IC2', 'IC3', 'IC4');
    colors <- annot_colors$icluster;
    foo <- foo[foo$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4'),]
  }

  plist[[i]] = ggboxplot(foo, category, 'score', xlab='', ylab='Dependency Score', title='ENRS', outlier.shape=NA,
                         facet.by='varlabel', scales='free_y', ncol=length(scores),
                         fill=category, palette=colors,  order=orders) +
    stat_compare_means(comparison=comps, size=5, label='p.signif') +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    theme(text = element_text(size = 12, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2)
}
print(ggarrange(plotlist=plist, nrow=2, ncol=1));
 
#-------------------- 
# (2) IC-PAM50 
#
# append data for multiple variables
for (i in 1:length(scores)){
  score <- scores[i];
  df <- data.frame(mdata, score=mdata[,score], varlabel=gsub('ENRS_CERES_', '', score));
  if (i == 1) { foo <- df; } else { foo <- rbind(foo, df)}
}

# reorder levels
foo$varlabel = factor(foo$varlabel, levels=c('ESR1', 'CDK4', 'CDK6', 'CDK2'))  

comps <- list(c('IC1:Basal', 'IC1:Her2'), c('IC1:Her2', 'IC1:LumB'), c('IC1:Basal', 'IC1:LumB'),
              c('IC1:LumB', 'IC2:LumB'), c('IC1:LumB', 'IC3:LumB'));
orders <- c('IC1:Basal',  'IC1:Her2', 'IC1:LumB', 'IC2:LumB', 'IC3:LumB');
colors <- annot_colors$ic_pam50;
foo <- foo[foo$ic_pam50 %in% c('IC1:Basal', 'IC1:Her2', 'IC1:LumB', 'IC2:LumB', 'IC3:LumB'),]

category = 'ic_pam50';
p = ggboxplot(foo, category, 'score', xlab='', ylab='Dependency Score', title='ENRS', outlier.shape=NA,
                       facet.by='varlabel', scales='free_y', ncol=length(scores), 
                       fill=category, palette=colors,  order=orders) +
  stat_compare_means(comparison=comps, size=5, label='p.signif') +
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  theme(text = element_text(size = 12, face='bold')) +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2) +
  rotate_x_text(45)

print(ggarrange(p, nrow=2, ncol=1, heights=c(2,1)));

#---------------------------------------------------------------------------
# (3) Longitudinal variable changes Pre vs. Post (paired BL/PT)
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

#---------------------------------------------------------
# (4) compare variable vs. key categories 
# gene MUT, AMP
#
genomic_muts = c('MUT.ESR1', 'MUT.RB1', 'MUT.TP53');
colors_alter = c('deepskyblue3', 'firebrick3');

mm <- melt(mdata[,c('sample_id', genomic_muts)], id.vars=1);
mm$value <- ifelse(mm$value == 1, 'MUT', 'WT');
mm$variable <- gsub('MUT.', '', mm$variable);

plist <- list();
comps <- list(c('WT', 'MUT'));
for (i in 1:length(scores)){
  yvar <- scores[i];
  foo <- merge(mm, mdata[,c('sample_id', yvar)], by='sample_id')
  plist[[i]] <- ggboxplot(foo, 'value', yvar, xlab='', ylab='Dependency Score', title=gsub('ENRS_CERES_', '', yvar), outlier.shape=NA, remove=NA,
                          fill='value', palette=colors_alter, facet.by='variable') + 
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    stat_compare_means(size=4, label='p.signif', comp=comps) + 
    theme(text = element_text(size = 12, face='bold')) + 
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2);
}

ggarrange(plotlist=plist, nrow=2, ncol=2, legend=F);


dev.off();
