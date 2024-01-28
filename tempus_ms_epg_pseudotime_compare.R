#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# June 2, 2023
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
outfile <- paste(outdir, 'tempus_ms_epg_pseudotime_integrated.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ms_epg_pseudotime_compare.pdf', sep='/');

# load Tempus features
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

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
comps = list(c('A', 'B'), c('A', 'C'), c('A', 'D'), c('C', 'D'));
ordered_branch = c('A', 'B', 'C', 'D', 'E')
colors_epg20 <- brewer.pal(5, 'Set1')[c(4,3,2,1,5)];
names(colors_epg20) <- ordered_branch;

# trend p-value using mann-kendall test
tmp <- tapply(mm$pseudotime, mm$branch_epg20, mean);
res <- mk.test(tmp[c('A', 'B', 'C', 'D')]);
title <- paste('EPG Branches (trend: p=', round(res$p.value, 5), ')', sep='');

p1 <- ggboxplot(mm, "branch_epg20", "pseudotime", xlab='EPG Branch', ylab='Pseudotime', 
                title=title, outlier.shape=NA,
                fill='branch_epg20', palette=colors_epg20,
                order=ordered_branch) +
  stat_compare_means(comparison=comps, size=4, label='p.format') +
  theme(text = element_text(size = 10, face='bold'), legend.position='none') +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

print(ggarrange(p1, nrow=2, ncol=2));

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
  tj_nodes <- as.character(tj_nodes);
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
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
}

print(ggarrange(plotlist=plist, nrow=2, ncol=2));

#-------------------------------------------------
# (4) association of EPG nodes vs. features
#
features <- c('EXPR.ESR1', 'EXPR.PGR', 'EXPR.CCNE1', 'EXPR.CCNE2', 'EXPR.CCND1', 'EXPR.CDK6', 'EXPR.SPDEF',
              'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response',
              'c2_cgp.SMID_BREAST_CANCER_ERBB2_UP', 'h.HALLMARK_ESTROGEN_RESPONSE_EARLY',
              'ENRS_CERES_ESR1', 'ENRS_CERES_CDK4', 'ENRS_CERES_CDK2')

c <- 1;
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(trajectories)){
  tj <- trajectories[i]; 
  tj_branches <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$branch;
  tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj,]$node;
  tj_nodes <- as.character(tj_nodes);
  colors_tj <- colors_epg20[tj_branches];
  names(colors_tj) <- tj_nodes;
  foo <- mm[mm$node_epg20 %in% tj_nodes,];

  for (j in 1:length(features)){
    feature <- features[j];
    ymin <- min(mm[,feature]);
    ymax <- max(mm[,feature]);
    
    # trend p-value using mann-kendall test
    tmp <- tapply(foo[, feature], foo$node_epg20, mean);
    res <- mk.test(tmp[tj_nodes])
    title <- paste(paste('Trajectory', tj, sep=': '), 
                   ' (trend: p=', round(res$p.value, 5), ')', sep='');
    
    plist[[c]] <- ggboxplot(foo, "node_epg20", feature, xlab='EPG Nodes', ylab=feature, 
                            ylim=c(ymin, ymax),
                          title=title, outlier.shape=NA,
                          fill='node_epg20', palette=colors_tj,
                          order=tj_nodes) +
    theme(text = element_text(size = 10, face='bold'), legend.position='none') +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
    c <- c + 1;
  }
}

print(ggarrange(plotlist=plist, nrow=2, ncol=2));

dev.off();

#---------------------------------------------------------------
# (5) heatmap of CSE drug resistance signatures (selected)
# sorted by pseudotime, EPG branch, PAM50 subtypes
#
library('pheatmap')
library('gplots')

rdata = tempus_bc();
sdata = rdata$sdata;

# load CSE GSVA data
gsvaobj <- tempus_bc_cse_gsva(rdata);
gsvadata <- gsvaobj$scores;
gscollections = gsvaobj$annotation$collection;
genesets <- c('h.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'h.HALLMARK_IL6_JAK_STAT3_SIGNALING',
              'h.HALLMARK_INFLAMMATORY_RESPONSE', 'h.HALLMARK_GLYCOLYSIS',
              'h.HALLMARK_MTORC1_SIGNALING', 'h.HALLMARK_CHOLESTEROL_HOMEOSTASIS',
              'h.HALLMARK_ESTROGEN_RESPONSE_LATE',
              'h.HALLMARK_ESTROGEN_RESPONSE_EARLY', 'h.MYC_TARGETS_V1', 'h.E2F_TARGETS',
              'h.HALLMARK_ALLOGRAFT_REJECTION', 'h.HALLMARK_INTERFERON_GAMMA_RESPONSE', 
              'h.HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'h.HALLMARK_HYPOXIA',
              'h.HALLMARK_IL2_STAT5_SIGNALING', 'h.HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY');
gsva <- gsvadata[rownames(gsvadata) %in% genesets,]
# match samples
gsva <- gsva[, mm$sample_id];

# sort signatures by pseudotime
sidx <- order(mm$pseudotime, decreasing=T);
foo <- gsva[,sidx];
mdata <- mm[sidx,];
rownames(foo) <- gsub('h.HALLMARK_', '', rownames(foo))

# column annotations
annot_cols <- data.frame(
  'pam50' = as.factor(mdata$subtype_pam50),
  'icluster' = as.factor(mdata$icluster),
  'prepost' = as.factor(mdata$pre_post),
  'pseudotime' = mdata$pseudotime
);
rownames(annot_cols) = mdata$sample_id;

# column annot. colors
annot_colors = list(
  prepost = c('Pre'='white', 'Post'='dimgrey'),
  pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  icluster = c('IC1'='forestgreen', 'IC2'='pink', 'IC3'='violetred', 'IC4'='purple4', 'IC5'='beige')
)

# under construction
pheatmap(foo, cex=1, cluster_cols=F, 
         annotation_col=annot_cols, annotation_colors=annot_colors, fontsize_row = 10,
         show_colnames=F, color=bluered(256), cellwidth=3, cellheight=10, border_color=NA,
         filename=paste(outdir, 'tempus_ms_peudotime_gsva_pheatmap.pdf', sep='/'));

#-------------------------------------
# key marker gene expression
#
sidx <- order(mm$pseudotime, decreasing=T);
features <- c('EXPR.PGR', 'EXPR.CCNE1', 'EXPR.ESR1');
foo <- t(mm[,features]);
colnames(foo) <- mm$sample_id;
foo <- foo[,sidx];

pheatmap(foo, cex=1, cluster_cols=F,
         annotation_col=annot_cols, annotation_colors=annot_colors, fontsize_row = 10,
         show_colnames=F, color=bluered(256), cellwidth=3, cellheight=10, border_color=NA,
         filename=paste(outdir, 'tempus_ms_peudotime_expr_pheatmap.pdf', sep='/'));



