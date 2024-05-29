#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 26, 2024
#
# Correlations of EPG trajectory vs. gene signatures
# (1) p-value, direction for trajectories A-B vs. A-C
# (2) p-value & correlation vs. pseudotime
# (3) identify signatures increased in A-C vs. A-B
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
outfile <- paste(outdir, 'tempus_ms_epg_trajectory_compare_gsva_cse.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ms_epg_trajectory_compare_gsva_cse.pdf', sep='/');

# load Tempus features
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

#--------------------------------
# Tempus BC R data object
# add gene signatures (CSE GSVA)
#
#rdata = tempus_bc();

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

# EPG tree with 3 branches & 20 nodes
comps = list(c('A', 'B'), c('B', 'C'), c('A', 'C'));
ordered_branch = c('A', 'B', 'C')
colors_epg20 <- brewer.pal(4, 'Set1')[c(4,3,2)];
names(colors_epg20) <- ordered_branch;

#---------------------------------------------------------------
# (2) calculate association of EPG nodes vs. signatures
#

# Hallmark genesets + GO BP "stem cell" related genests
idx <- (rdata$Analysis$GSVA_annotation$pathway_class == 'hallmark gene sets');
#idx <- (rdata$Analysis$GSVA_annotation$pathway_class == 'hallmark gene sets') |
#  ((rdata$Analysis$GSVA_annotation$pathway_class == 'C5 BP: GO biological process' & 
#      grepl('STEM_CELL', rdata$Analysis$GSVA_annotation$pathway_name)));
genesets <- rdata$Analysis$GSVA_annotation$pathway_name[idx];

# extract CSE GSVA values
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

# correlation vs. pseudotime
#
df1 <- NULL;
for (i in 1:length(genesets)){
  geneset <- genesets[i];
  res <- cor.test(mm$pseudotime, mm[,geneset]);
  df1 <- rbind(df1, data.frame(geneset, correlation.PT=res$estimate, pval.PT=res$p.value))
}

# trend p-value using mann-kendall test
# direction - increase vs. decrease
#
graph <- 'EPG20';
df <- NULL;
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph])
for (i in 1:length(genesets)){
  
  geneset <- genesets[i];
  tj1 <- trajectories[1];
  tj2 <- trajectories[2];

  # trajectory A-B
  tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj1,]$node;
  tj_nodes <- as.character(unique(tj_nodes));
  foo <- mm[mm$node_epg20 %in% tj_nodes,];
  tmp <- tapply(foo[, geneset], foo$node_epg20, mean, na.rm=T);
  res1 <- mk.test(tmp[tj_nodes]);

  tj_nodes <- epgdata_trajectory[epgdata_trajectory$graph == graph & epgdata_trajectory$trajectory == tj2,]$node;
  tj_nodes <- as.character(unique(tj_nodes));
  foo <- mm[mm$node_epg20 %in% tj_nodes,];
  tmp <- tapply(foo[, geneset], foo$node_epg20, mean, na.rm=T);
  res2 <- mk.test(tmp[tj_nodes]);
  
  dd <- data.frame(geneset, zval=res2$statistic, pval=res2$p.value, 
                   logp=ifelse(res2$statistic>0, -log2(res2$p.value), log2(res2$p.value)),
                   zval2=res1$statistic, pval2=res1$p.value, 
                   logp2=ifelse(res1$statistic>0, -log2(res1$p.value), log2(res1$p.value)));
  df <- rbind(df, dd);
}
df$lods <- df$logp - df$logp2;
df$qval <- p.adjust(df$pval, method='BH')


# geneset classification
gs_metabolism <- c('GLYCOLYSIS', 'CHOLESTEROL_HOMEOSTASIS', 'MTORC1_SIGNALING');
gs_growth_proliferation <- c('G2M_CHECKPOINT', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'MYC_TARGETS_V2');
gs_inflammation <- c('IL6_JAK_STAT3_SIGNALING', 'IL2_STAT5_SIGNALING', 
                     'INFLAMMATORY_RESPONSE', 'ALLOGRAFT_REJECTION', 'TNFA_SIGNALING_VIA_NFKB');
#gs_stemcell <- c('APICAL_JUNCTION')


#-------------------------------------------
# (3) create plots for significant hits
# compare trajectories side-by-side
#

# graphics output file
pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);

#genesets <- df$geneset[df$lods > 10 & df$qval < 0.05];
genesets <- df$geneset[df$lods > 10];

c <- 1;
plist <- list();
graph <- 'EPG20';
colors_branch <- rainbow(5);
trajectories <- unique(epgdata_trajectory$trajectory[epgdata_trajectory$graph == graph]);
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


#-------------------------------------------------
# (4) plot trajectory for selected signatures
#

genesets <- c(
  'c2_cgp.SMID_BREAST_CANCER_ERBB2_UP',   'h.HALLMARK_ESTROGEN_RESPONSE_EARLY',
  'h.HALLMARK_G2M_CHECKPOINT', 'h.HALLMARK_E2F_TARGETS', 'h.HALLMARK_MYC_TARGETS_V1',
    'c2cp.WP_HIPPOYAP_SIGNALING_PATHWAY',
  'c5_bp.GO_SOMATIC_STEM_CELL_DIVISION',  'c5_bp.GO_STEM_CELL_DIFFERENTIATION',
  'c5_bp.GO_REGULATION_OF_STEM_CELL_DIFFERENTIATION',
  'h.HALLMARK_NOTCH_SIGNALING', 'h.HALLMARK_CHOLESTEROL_HOMEOSTASIS',
  'h.HALLMARK_APICAL_JUNCTION', 'h.HALLMARK_INFLAMMATORY_RESPONSE', 'h.HALLMARK_GLYCOLYSIS',
  'h.HALLMARK_IL6_JAK_STAT3_SIGNALING', 'h.HALLMARK_KRAS_SIGNALING_UP')
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


#-----------------------------------------------------------------
# (5) Compare gene signature "log-odds", Trajector A-C vs. A-B
# scatter plot, colored by DE status
#
library('RColorBrewer');
colors_epg_traj <- brewer.pal(5, 'Set1')[c(3,1,2,4)];
colors_epg_traj[3] <- colors_epg_traj[4] <- 'grey';

# Graphics output file
#plotfile2 <- paste(outdir, paste('tempus_ms_epg_trajectory_gsva_scatter.pdf', sep='.'), sep='/');
#pdf(file=plotfile2, onefile=T, height=10, width=12, pointsize=12);

foo <- df;
logpcut <- 5;
logpcut2 <- 3;
foo$traj_enrichment <- ifelse((foo$logp > logpcut & foo$logp2 > logpcut), 'Both',
                                  ifelse(foo$logp > logpcut & foo$logp2 < logpcut2, 'A-C',
                                         ifelse(foo$logp2 > logpcut & foo$logp2 < logpcut2, 'A-B', 'Neither')));

selected <- c('INFLAMMATORY_RESPONSE', 'TNFA_SIGNALING_VIA_NFKB', 'APICAL_JUNCTION',
              'APICAL_JUNCTION', 'P53_PATHWAY',
              'IL6_JAK_STAT3_SIGNALING', 'IL2_STAT5_SIGNALING', 'KRAS_SIGNALING_UP',
              'ESTROGEN_RESPONSE_EARLY', 'EPITHELIAL_MESENCHYMAL_TRANSITION');

# scatter-plot v2
ggscatter(foo, x='logp', y='logp2', 
          xlim=c(-12, 15), ylim=c(-12, 15),
          xlab='Log2p (Trajectory A-C Trend)', ylab='Log2p (Trajectory A-B Trend)', 
          color='traj_enrichment', palette=colors_epg_traj,
          size='logp', alpha=.7,
          label='geneset', repel=TRUE,
          label.select=selected) +
  #stat_cor(method='spearman', label.x=-5, label.y=15, size=7) +
  scale_size(range = c(.1, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  theme(text = element_text(size = 20, face='bold'), legend.position='right') 

dev.off()

#----------------------------------------------------
# output - PT, EPG branch, node, trajectories
#
out <- merge(foo, df1, by='geneset', sort=F, all.x=T)
write.table(out, file=outfile, sep="\t", col.names=T, row.names=F, quote=F);


