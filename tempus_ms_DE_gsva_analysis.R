#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# July 03, 2023
#
# script to integrate hits from multiple DE comparisons
# comparisons: PAM50, Pre/Post, IC1-4, PFS assoc
# data type: CSE, bulk expression
# variable type: gene expression, gene signature (gsva)
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library('pheatmap');
library('gplots');
library(EnhancedVolcano)

file1 <- paste(datadir, 'DE/tempus_cse_DE_GSVAs_pam50_wo_PI.txt', sep='/');
file2 <- paste(datadir, 'DE/tempus_cse_DE_GSVAs_iclusters_wo_PI.txt', sep='/');
file3 <- paste(datadir, 'DE/tempus_cse_DE_GSVAs_prepost_wo_PI.txt', sep='/');
file4 <- paste(datadir, 'DE/tempus_cse_pfs_association_GSVA.txt', sep='/');
  
outfile <- paste(outdir, 'tempus_ms_DE_gsva_analysis.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ms_DE_gsva_analysis.pdf', sep='/');

#------------------------------------------------------------------
# (1) parse in PAM50 DE analysis results (multinomial)
#
d1 <- read.table(file1, header=T, sep="\t");
colnames(d1) <- tolower(colnames(d1));
groups <- c('basal', 'her2', 'luma', 'lumb');
for (i in 1:length(groups)){
  group <- groups[i];
  pval <- d1[,paste(group, 'pval', sep='_')];
  logp <- ifelse(d1[,group] > 0, -log10(pval), log10(pval));
  d1[, paste(group, 'logp', sep='_')] <- logp;
}
idx <- colnames(d1) %in% groups;
colnames(d1)[idx] <- paste(colnames(d1)[idx], 'coef', sep='_');

#------------------------------------------------------------------
# (2) parse in IC1-4 DE analysis results (multi-nomial)
#
d2 <- read.table(file2, header=T, sep="\t");
colnames(d2) <- tolower(colnames(d2));
groups <- c('ic1', 'ic2', 'ic3', 'ic4');
for (i in 1:length(groups)){
  group <- groups[i];
  pval <- d2[,paste(group, 'pval', sep='_')];
  logp <- ifelse(d2[,group] > 0, -log10(pval), log10(pval));
  d2[, paste(group, 'logp', sep='_')] <- logp;
}
idx <- colnames(d2) %in% groups;
colnames(d2)[idx] <- paste(colnames(d2)[idx], 'coef', sep='_');

#------------------------------------------------------------------
# (3) parse in Pre/Post DE analysis results (binomial)
#
d3 <- read.table(file3, header=T, sep="\t");
colnames(d3)[-1] <- c('prepost_logfc', 'prepost_coef', 'prepost_fval', 'prepost_pval', 'prepost_logp', 'prepost_qval')

#------------------------------------------------------------------
# (4) parse in PFS assoc analysis results (binomial)
#
d4 <- read.table(file4, header=T, sep="\t");
d4$logp <- ifelse(d4$hr>1, -log10(d4$pval), log10(d4$pval))
d4m <- d4[,c('geneset', 'statistic', 'hr', 'pval', 'qval', 'logp')];

#------------------------------------------------------------------
# (5) integrate across comparisons into table
#
d12 <- merge(d1, d2, by='geneset', all.x=T);
d123 <- merge(d12, d3, by='geneset', all.x=T);
d1234 <- merge(d123, d4m, by='geneset', all.x=T);

# rearrange columns
hdrs <- c('geneset', 
          c(colnames(d1234)[grepl('lods', colnames(d1234))],
            colnames(d1234)[grepl('logp', colnames(d1234))],
            colnames(d1234)[grepl('coef', colnames(d1234))],
            colnames(d1234)[grepl('pval', colnames(d1234))],
            colnames(d1234)[grepl('qval', colnames(d1234))],
            colnames(d1234)[grepl('logfc', colnames(d1234))],
            colnames(d1234)[grepl('hr', colnames(d1234))]));
mdata <- d1234[grepl('^h.HALLMARK', d1234$geneset), hdrs];
mdata$geneset <- gsub('h.HALLMARK_', '', mdata$geneset);

#----------------------------------------
# (6) Volcano plot
#
# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=14, pointsize=12);

foo <- d1[grepl('^h.HALLMARK', d1$geneset),]
foo$geneset <- gsub('h.HALLMARK_', '', foo$geneset)

# Basal DE genesets
EnhancedVolcano(foo, lab=foo$geneset, pointSize = 8, title='Basal',
                x='basal_coef', y='basal_pval',
                FCcutoff=.3, pCutoff=0.05,
                xlim=c(-1, 1), ylim=c(0,4),
                selectLab = c('IL6_JAK_STAT3_SIGNALING', 'INFLMMATORY_RESPONSE',
                              'MYC_TARGETS_V2', 'MYC_TARGETS_V1', 'GLYCOLYSIS',
                              'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INTERFERON_GAMMA_RESPONSE',
                              'ALLOGRAFT_REJECTION', 'ESTROGEN_RESPONSE_EARLY'),
                axisLabSize = 20, labSize=7, legendLabSize = 20,
                boxedLabels = T, labFace='bold', xlab = 'Coef',
                drawConnectors=TRUE, colConnectors='black', widthConnectors=1,
                gridlines.major = F, gridlines.minor = F,
                legendLabels = c("NS", "Coef", "p-value", expression("p-value" ~ and ~ "Coef")),)

# Her2 DE genesets
EnhancedVolcano(foo, lab=foo$geneset, pointSize = 8, title='Her2',
                x='her2_coef', y='her2_pval',
                FCcutoff=.3, pCutoff=0.05,
                xlim=c(-1, 1), ylim=c(0,3),
                selectLab = c('PI3K_AKT_MTOR_SIGNALING', 'MTORC1_SIGNALING', 'GLYCOLYSIS', 
                              'OXIDATIVE_PHOSPHORYLATION', 'ESTROGEN_RESPONSE_EARLY'),
                axisLabSize = 20, labSize=7, legendLabSize = 20,
                boxedLabels = T, labFace='bold', xlab = 'Coef',
                drawConnectors=TRUE, colConnectors='black', widthConnectors=1,
                gridlines.major = F, gridlines.minor = F,
                legendLabels = c("NS", "Coef", "p-value", expression("p-value" ~ and ~ "Coef")),)

foo <- d2[grepl('^h.HALLMARK', d2$geneset),]
foo$geneset <- gsub('h.HALLMARK_', '', foo$geneset)

# IC1 DE genesets
EnhancedVolcano(foo, lab=foo$geneset, pointSize = 8, title='IC1',
                x='ic1_coef', y='ic1_pval',
                FCcutoff=.3, pCutoff=0.05,
                xlim=c(-1, 1), ylim=c(0,8),
                selectLab = c('E2F_TARGETS', 'G2M_CHECKPOINT', 'EPITHELIAL_MESENCHYMAL_TRANSITION',
                              'MTORC1_SIGNALING', 'MYC_TARGETS_V2', 'MYC_TARGETS_V1',
                              'GLYCOLYSIS', 'ESTROGEN_RESPONSE_EARLY'),
                axisLabSize = 20, labSize=7, legendLabSize = 20,
                boxedLabels = T, labFace='bold', xlab = 'Coef',
                drawConnectors=TRUE, colConnectors='black', widthConnectors=1,
                gridlines.major = F, gridlines.minor = F,
                legendLabels = c("NS", "Coef", "p-value", expression("p-value" ~ and ~ "Coef")),)


#----------------------------------------
# (6) heatmap with waterfall side
#
library('ComplexHeatmap');

rankedby <- c('ic1_logp', 'basal_logp', 'her2_logp');
nrows = 10;
for (i in 1:3){
  rb <- rankedby[i];
  idx1 <- order(mdata[, rb], decreasing=T)[1:nrows];
  idx2 <- order(mdata[, rb], decreasing=F)[1:nrows];

  colnames <- colnames(mdata)[grepl('logp', colnames(mdata))][c(5:8, 1:4, 9:10)]
  foo <- mdata[c(idx1, rev(idx2)), colnames];
  rownames(foo) <- mdata$geneset[c(idx1, rev(idx2))];
  colnames(foo) <- c('IC1', 'IC2', 'IC3', 'IC4', 'Basal', 'Her2', 'LumA', 'LumB', 'PrePost', 'PFS_Assoc')
  foo <- foo[,c(1:2,4,3,5:10)]
  top_ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 'white'),
                              labels = c("ICLUSTER", "PAM50"),
                              labels_gp = gpar(col='black', fontsize=14, fontface='bold')));
  
  # Pre/Post p-value
  barcols <- c('blue', 'red')[as.numeric(as.factor(foo$PrePost > 0))]
  row_ha <- rowAnnotation(PrePost=anno_barplot(foo$PrePost, gp=gpar(fill=barcols)),
                          width=unit(30,'mm'));
  
  # PFS assoc p-value
  barcols2 <- c('blue', 'red')[as.numeric(as.factor(foo$PFS_Assoc > 0))]
  row_ha2 <- rowAnnotation(PFS_Assoc=anno_barplot(foo$PFS_Assoc, gp=gpar(fill=barcols2)),
                          width=unit(30,'mm'));
  p <- Heatmap(as.matrix(foo[,1:8]), cluster_columns=F, cluster_rows=F, 
          name='-log10(p)', 
          right_annotation = c(row_ha, row_ha2),
          column_split=rep(c('1', '2'), each=4), column_gap=unit(5,'mm'),
          row_split=rep(c('1', '2'), each=nrows), row_gap=unit(3, 'mm'),
          column_title=NULL, row_title=NULL,
          top_annotation = top_ha,
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 15),
          width = ncol(foo[,1:8])*unit(9, "mm"),
          height = nrow(foo[,1:8])*unit(9, "mm"))
  
  print(p);
}

dev.off();






  