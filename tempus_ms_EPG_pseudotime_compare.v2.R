#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 21, 2024
#
# Quick characterization analysis of EPG branch, node and trajectory
# EPG20 (20-node tree); EPG30 (30-node tree)
#
source('Code/TEMPUS_BC.R');
library('ggpubr');
library('RColorBrewer');
library('trend');

datadir = 'Input';
outdir = 'Output.Manuscript';

# input files
infile <- paste(datadir, '8.sample.clustering.output.batchCorrected.wPseudotime.csv', sep='/');
infile2 <- paste(datadir, 'tempus_epg_node_branch.txt', sep='/');
infile3 <- paste(datadir, 'tempus_epg_trajectory.txt', sep='/');
outfile <- paste(outdir, 'tempus_ms_epg_pseudotime_integrated.v2.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ms_EPG_pseudotime_compare.v2.pdf', sep='/');
plotfile2 <- paste(outdir, 'tempus_ms_epg_trajectory_compare_gsva_cse.v2.pdf', sep='/');

# load Tempus features
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

#--------------------------------
# Tempus BC R data object
# add gene signatures (CSE GSVA)
#
# rdata = tempus_bc();

#--------------------------------------
# (1) load Monocle "pseudotime"
#
indata <- read.csv(infile);
ptdata <- indata[,c('X', 'pseudotime')];
colnames(ptdata) <- c('sample_id', 'pseudotime');

# load EPG branch & node data
epgdata <- read.table(infile2, header=T, sep='\t');

# load EPG trajectory data
epgdata_trajectory <- read.table(infile3, header=T, sep='\t');

# merge data
m1 <- merge(ptdata, epgdata, by='sample_id')
mm <- merge(m1, mdata, by='sample_id')

mm$branch_epg20 <- mm$branch;
mm$node_epg20 <- mm$node;

#----------------------------------------------------
# output - PT, EPG branch, node, trajectories
#
out <- subset(mm, select=c('sample_id', 'patient_id', 'pseudotime', 'branch_epg20', 'node_epg20'))
write.table(out, file=outfile, sep="\t", col.names=T, row.names=F, quote=F);

#--------------------------------------------
# (2) association of EPG branch vs. PT
#

# graphics output file
pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);

# EPH tree with 20 nodes, 5 branches
comps = list(c('A', 'B'), c('B', 'C'), c('A', 'C'));
ordered_branch = c('A', 'B', 'C')
colors_epg20 <- brewer.pal(4, 'Set1')[c(4,3,2)];
names(colors_epg20) <- ordered_branch;


# increase %IC1 vs. EPG nodes along the trajectory
colors_icluster = brewer.pal(5, 'Set1')
foo <- mm[mm$icluster != 'IC5',];
par(mfrow=c(2,1));
barplot(prop.table(table(foo$icluster, foo$node_epg20), 2), col=colors_icluster, ylim=c(0, 1));
frame();
legend('top', fill=colors_icluster, legend=levels(as.factor(mm$icluster)))


# trend p-value using mann-kendall test
tmp <- tapply(mm$pseudotime, mm$branch_epg20, mean);
res <- mk.test(tmp[c('A', 'B', 'C')]);
title <- paste('EPG Branches (trend: p=', round(res$p.value, 5), ')', sep='');

p1 <- ggboxplot(mm, "branch_epg20", "pseudotime", xlab='EPG Branch', ylab='Pseudotime', 
                title=title, outlier.shape=NA,
                fill='branch_epg20', palette=colors_epg20,
                order=ordered_branch) +
  stat_compare_means(comparison=comps, size=6, label='p.signif') +
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  theme(text = element_text(size = 10, face='bold'), legend.position='none') +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

print(ggarrange(p1, nrow=2, ncol=3));

#--------------------------------------
# (3) association of EPH nodes vs. PT
#

# EPG20 trajectories
#
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
  foo <- mm[mm$node_epg20 %in% tj_nodes,];
  
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
features <- c('EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1', 'EXPR.CCNE2', 'EXPR.CCND1', 'EXPR.CDK6', 'EXPR.SPDEF', 'EXPR.FOXA1',
              'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response',
              'c2_cgp.SMID_BREAST_CANCER_ERBB2_UP', 'h.HALLMARK_ESTROGEN_RESPONSE_EARLY',
              'ENRS_CERES_ESR1', 'ENRS_CERES_CDK4', 'ENRS_CERES_CDK2', 'mutsig_13_sigma')

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
  foo <- mm[mm$node_epg20 %in% tj_nodes,];
  
  for (j in 1:length(features)){
    feature <- features[j];
    ymin <- min(mm[,feature], na.rm=T);
    ymax <- max(mm[,feature], na.mr=T);
    
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

#-------------------------------------------------
# (5) association of EPG nodes vs. signatures
#

# graphics output file
pdf(file=plotfile2, onefile=T, height=8, width=10, pointsize=12);

genesets <- c(
    'c2_cgp.SMID_BREAST_CANCER_ERBB2_UP', 'h.HALLMARK_ESTROGEN_RESPONSE_EARLY',
    'c2cp.WP_HIPPOYAP_SIGNALING_PATHWAY',
    'c5_bp.GO_SOMATIC_STEM_CELL_DIVISION',
    'c5_bp.GO_REGULATION_OF_STEM_CELL_DIFFERENTIATION',
    'h.HALLMARK_NOTCH_SIGNALING', 'h.HALLMARK_CHOLESTEROL_HOMEOSTASIS',
    'h.HALLMARK_APICAL_JUNCTION', 'h.HALLMARK_INFLAMMATORY_RESPONSE', 'h.HALLMARK_GLYCOLYSIS',
    'h.HALLMARK_IL6_JAK_STAT3_SIGNALING', 
    'h.HALLMARK_KRAS_SIGNALING_UP', 'h.HALLMARK_MYC_TARGETS_V1')
tmp <- rdata$Analysis$GSVA_BayesPrism_cancer[genesets, mdata$sample_id]

rownames(tmp) <- gsub('c2_cgp.', '', rownames(tmp));
rownames(tmp) <- gsub('c5_bp.', '', rownames(tmp));
rownames(tmp) <- gsub('c2cp.', '', rownames(tmp));
rownames(tmp) <- gsub('h.HALLMARK_', '', rownames(tmp));
genesets <- rownames(tmp);
mdata <- cbind(mdata, t(tmp))
mm <- merge(m1, mdata, by='sample_id')
mm$branch_epg20 <- mm$branch;
mm$node_epg20 <- mm$node;

c <- 1;
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(genesets)){
  geneset <- genesets[i];
  ymin <- min(mm[,geneset], na.rm=T);
  ymax <- max(mm[,geneset], na.rm=T);
  
  for (j in 1:length(trajectories)){
    tj <- trajectories[j]; 
    tj_branches <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$branch;
    tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$node;
    tj_nodes <- as.character(unique(tj_nodes));
    colors_tj <- colors_epg20[tj_branches];
    names(colors_tj) <- tj_nodes;
    foo <- mm[mm$node_epg20 %in% tj_nodes,];
    
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



