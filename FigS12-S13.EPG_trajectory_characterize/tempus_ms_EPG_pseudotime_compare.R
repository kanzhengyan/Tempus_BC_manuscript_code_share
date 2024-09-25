#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 21, 2024
# Sept. 21, 2024
#
# Quick characterization analysis of EPG branch, node and trajectory
#
setwd('../FigS12-S13.EPG_trajectory_characterize');
library('ggpubr');
library('RColorBrewer');
library('trend');

datadir = 'Input';
outdir = 'Output.Manuscript';

# input files
infile <- paste(datadir, 'tempus_epg_trajectory.txt', sep='/');
plotfile <- paste(outdir, 'Fig-S12.EPG_pseudotime_compare.pdf', sep='/');
plotfile2 <- paste(outdir, 'Fig-S13.EPG_trajectory_compare_gsva_cse.pdf', sep='/');

# load Tempus molecular features
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

# Monocle "pseudotime"
ptdata <- mdata[, c('sample_id', 'pseudotime')]; 

# EPG trajectory data
epgdata_trajectory <- read.table(infile, header=T, sep='\t');

# load CSE GSVA data
gsvafile <- paste(datadir, 'tempus_ms_supp_data.hallmark_cse_gsva.txt', sep='/');
gsvadata <- read.table(gsvafile, sep='\t', header=T)
rownames(gsvadata) <- gsvadata[,1]
gsva <- t(gsvadata[,-1]);

#-------------------------------------------------------------------
# (1) association of EPG trajectory vs. disease progression
#

# graphics output file
pdf(file=plotfile, onefile=T, height=8, width=12, pointsize=12);

# EPH tree with 20 nodes, 5 branches
comps = list(c('A', 'B'), c('B', 'C'), c('A', 'C'));
ordered_branch = c('A', 'B', 'C')
colors_epg20 <- brewer.pal(4, 'Set1')[c(4,3,2)];
names(colors_epg20) <- ordered_branch;
colors_icluster = c(brewer.pal(4, 'Set1'), 'grey');
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');
colors_prepost_pfe = c('Pre'='lightgrey', 'During'='grey60', 'Post'='grey20');

foo <- mdata;
foo$pre_post <- factor(foo$pre_post, levels=c('Pre', 'Post'))
foo$pre_post_pfe <- factor(foo$pre_post_pfe, levels=c('Pre', 'During', 'Post'))

foo1 <- foo[foo$branch_epg20 %in% c('A', 'B'), ]
foo2 <- foo[foo$branch_epg20 %in% c('A', 'C'), ]

par(mfrow=c(2,2));
barplot(prop.table(table(foo1$pre_post_pfe, foo1$node_epg20), 2), col=colors_prepost_pfe, 
        ylim=c(0, 1), ylab='% samples', main='Trajectory A-B');
barplot(prop.table(table(foo2$pre_post_pfe, foo2$node_epg20), 2), col=colors_prepost_pfe, 
        ylim=c(0, 1), ylab='% samples', main='Trajectory A-C');

barplot(prop.table(table(foo1$icluster, foo1$node_epg20), 2), col=colors_icluster, 
        ylim=c(0, 1), ylab='% samples', main='Trajectory A-B');
barplot(prop.table(table(foo2$icluster, foo2$node_epg20), 2), col=colors_icluster, 
        ylim=c(0, 1), ylab='% samples', main='Trajectory A-C');

frame();
legend('top', fill=colors_icluster, legend=levels(as.factor(mdata$icluster)), ncol=5)
legend('center', fill=colors_prepost_pfe, legend=levels(as.factor(foo$pre_post_pfe)), ncol=3);
legend('bottom', fill=colors_prepost, legend=levels(as.factor(foo$pre_post)), ncol=2)

#-------------------------------------------------
# (2) association of EPH branch/nodes vs. PT
# trend p-value by mann-kendall test
#
tmp <- tapply(mdata$pseudotime, mdata$branch_epg20, mean);
res <- mk.test(tmp[c('A', 'B', 'C')]);
title <- paste('EPG Branches (trend: p=', round(res$p.value, 5), ')', sep='');

p1 <- ggboxplot(mdata[mdata$branch_epg20 %in% ordered_branch,], "branch_epg20", "pseudotime", xlab='EPG Branch', ylab='Pseudotime', 
                title=title, outlier.shape=NA,
                fill='branch_epg20', palette=colors_epg20,
                order=ordered_branch) +
  stat_compare_means(comparison=comps, size=6, label='p.signif') +
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  theme(text = element_text(size = 10, face='bold'), legend.position='none') +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

print(ggarrange(p1, nrow=2, ncol=3));

# EPG20 trajectories
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(trajectories)){
  tj <- trajectories[i]; 
  tj_branches <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$branch;
  tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$node;
  tj_nodes <- as.character(unique(tj_nodes));
  colors_tj <- colors_epg20[tj_branches]
  names(colors_tj) <- tj_nodes;
  foo <- mdata[mdata$node_epg20 %in% tj_nodes,];
  
  # trend p-value using mann-kendall test
  tmp <- tapply(foo$pseudotime, foo$node_epg20, mean);
  res <- mk.test(tmp[tj_nodes])
  title <- paste(paste('Trajectory', tj, sep=': '), 
                 ' (trend: p=', round(res$p.value, 5), ')', sep='');
  
  plist[[i]] <- ggboxplot(foo, "node_epg20", "pseudotime", xlab='EPG Nodes', ylab='Pseudotime', 
                          title=title, outlier.shape=NA,
                          fill='node_epg20', palette=colors_tj,
                          order=tj_nodes) +
    theme(text = element_text(size = 10, face='bold'), legend.position='none') +
    geom_smooth(method='loess', se=T, aes(group=1, col='grey')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
}
print(ggarrange(plotlist=plist, nrow=2, ncol=2));

#-------------------------------------------------
# (4) association of EPG nodes vs. features
#
features <- c('EXPR.ESR1')

c <- 1;
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(trajectories)){
  tj <- trajectories[i]; 
  tj_branches <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$branch;
  tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$node;
  tj_nodes <- as.character(unique(tj_nodes));
  colors_tj <- colors_epg20[tj_branches];
  names(colors_tj) <- tj_nodes;
  foo <- mdata[mdata$node_epg20 %in% tj_nodes,];
  
  for (j in 1:length(features)){
    feature <- features[j];
    ymin <- min(mdata[,feature], na.rm=T);
    ymax <- max(mdata[,feature], na.mr=T);
    
    # trend p-value using mann-kendall test
    tmp <- tapply(foo[, feature], foo$node_epg20, mean, na.rm=T);
    res <- mk.test(tmp[tj_nodes])
    title <- paste(paste('Trajectory', tj, sep=': '), 
                   ' (trend: p=', round(res$p.value, 5), ')', sep='');
    
    plist[[c]] <- ggboxplot(foo, "node_epg20", feature, xlab='EPG Nodes', ylab=feature, 
                            ylim=c(ymin, ymax),
                            title=title, outlier.shape=NA,
                            fill='node_epg20', palette=colors_tj,
                            order=tj_nodes) +
      theme(text = element_text(size = 10, face='bold'), legend.position='none') +
      geom_smooth(method='loess', se=T, aes(group=1, col='grey')) +
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
    c <- c + 1;
  }
}

print(ggarrange(plotlist=plist, nrow=2, ncol=2));

dev.off();

#-----------------------------------------------------
# (5) association of EPG nodes vs. CSE signatures
#

# graphics output file
pdf(file=plotfile2, onefile=T, height=8, width=10, pointsize=12);

rownames(gsva) <- gsub('c2_cgp.', '', rownames(gsva));
rownames(gsva) <- gsub('c5_bp.', '', rownames(gsva));
rownames(gsva) <- gsub('c2cp.', '', rownames(gsva));
genesets <- c('MYC_TARGETS_V1', 'GLYCOLYSIS',
              'IL6_JAK_STAT3_SIGNALING', 'KRAS_SIGNALING_UP',
              'WP_HIPPOYAP_SIGNALING_PATHWAY',
              'GO_SOMATIC_STEM_CELL_DIVISION')
tmp <- gsva[genesets, mdata$sample_id]
mdata <- cbind(mdata, t(tmp));

c <- 1;
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(genesets)){
  geneset <- genesets[i];
  ymin <- min(mdata[,geneset], na.rm=T);
  ymax <- max(mdata[,geneset], na.rm=T);
  
  for (j in 1:length(trajectories)){
    tj <- trajectories[j]; 
    tj_branches <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$branch;
    tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$node;
    tj_nodes <- as.character(unique(tj_nodes));
    colors_tj <- colors_epg20[tj_branches];
    names(colors_tj) <- tj_nodes;
    foo <- mdata[mdata$node_epg20 %in% tj_nodes,];
    
    # trend p-value using mann-kendall test
    tmp <- tapply(foo[, geneset], foo$node_epg20, mean, na.rm=T);
    res <- mk.test(tmp[tj_nodes])
    title <- paste(paste('Trajectory', tj, sep=': '), 
                   ' (trend: p=', round(res$p.value, 5), ')', sep='');
    
    plist[[c]] <- ggboxplot(foo, "node_epg20", geneset, xlab='EPG Nodes', ylab=geneset, 
                            ylim=c(ymin, ymax),
                            title=title, outlier.shape=NA,
                            fill='node_epg20', palette=colors_tj,
                            order=tj_nodes) +
      theme(text = element_text(size=10, face='bold'), legend.position='none') +
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3) +
      geom_smooth(method='loess', se=T, aes(group=1, col='grey'))
    c <- c + 1;
  }
}
print(ggarrange(plotlist=plist, nrow=2, ncol=2));

dev.off();



