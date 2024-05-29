#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# Oct. 7, 2023
#
# show molecular covariates of Pseudotime & EPG branch as heatmap
# highlight CSE drug resistance signatures
# sorted by Pseudotime
#
source('Code/TEMPUS_BC.R');
library('ggpubr');
library('RColorBrewer');
library('ComplexHeatmap');

datadir = 'Input';
outdir = 'Output.Manuscript';

# input files
infile <- paste(datadir, 'tempus_ms_epg_pseudotime_integrated.txt', sep='/');
outfile <- paste(outdir, 'tempus_ms_pseudotime_covariate_analysis.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ms_pseudotime_covariate_heatmap.pdf', sep='/');

#--------------------------------------
# (1) load data 
# - "pseudotime" and "EPG" data
# - Hallmark CSE gene signatures
# - molecular features

# load Tempus data
rdata = tempus_bc();
sdata = rdata$sdata;

# load EPG branch & node data
indata <- read.table(infile, header=T, sep='\t');

# load Tempus features
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended_v3.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');

# merge feature vs. PT data
mm <- merge(mdata, indata, by='sample_id')

# load CSE GSVA data
gsvaobj <- tempus_bc_cse_gsva(rdata);
gsvadata <- gsvaobj$scores;
gscollections = gsvaobj$annotation$collection;

#-------------------------------------
# (2) select covariates for heatmap
#

# select genesets based on correlation vs. PT
genesets <- rownames(gsvadata)[grepl('h.HALLMARK', rownames(gsvadata))];
gsva <- gsvadata[rownames(gsvadata) %in% genesets, mm$sample_id];
cc <- cor(t(rbind(gsva, mm$pseudotime)));
ccvals <- cc[51,1:50]

pvals <- apply(gsva, 1, function(x) {
  res <- cor.test(x, mm$pseudotime);
  return(res$p.value)
})
rvals <- apply(gsva, 1, function(x) {
  res <- cor.test(x, mm$pseudotime, method='pearson');
  return(res$estimate)
})
res <- data.frame(geneset=names(pvals), p=pvals, r=rvals);
res$q <- p.adjust(res$p, method='BH');

genesets_cor <- res$geneset[abs(res$r) > 0.35]
gsva_cor <- gsva[rownames(gsva) %in% genesets_cor,];

write.table(res, file=outfile, sep="\t", col.names=T, row.names=F, quote=F);

# 2.1 Hallmark CSE gene signatures

# sort signatures by pseudotime
sidx <- order(mm$pseudotime, decreasing=T);
foo1 <- gsva_cor[,sidx];
rownames(foo1) <- gsub('h.HALLMARK_', '', rownames(foo1))

# 2.2 key marker gene expression
features <- c('EXPR.PGR', 'EXPR.CCNE1', 'EXPR.ESR1');
foo2 <- t(mm[sidx,features]);
colnames(foo2) <- mm[sidx,]$sample_id;

# 2.3 Molecular features (numerical)
features <- c('Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response',
              'Paloma3_F7_IFNG_response', 'Paloma3_F8_EMT');
foo3 <- t(scale(mm[sidx,features]));
colnames(foo3) <- mm[sidx,]$sample_id;
rownames(foo3) <- gsub('Paloma3_', '', rownames(foo3))

# 2.4 PAM50
features <- c('pam50_cor_luma', 'pam50_cor_lumb', 'pam50_cor_her2', 'pam50_cor_basal');
foo4 <- t(scale(mm[sidx,features]));
colnames(foo4) <- mm[sidx,]$sample_id;

#----------------------------------------
# (3) create complex heatmap
#
# columns sorted by Pseudotime (barplot)
# - key categories
# - CSE GSVA
# - marker gene expression
# - molecular features
#

# Graphics output file
pdf(file=plotfile, onefile=T, height=15, width=33, pointsize=12);

colors_set1 = brewer.pal(5, 'Set1')

# column annot. colors
annot_colors = list(
  prepost = c('Pre'='white', 'Post'='dimgrey'),
  pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  icluster = c('IC1'='forestgreen', 'IC2'='pink', 'IC3'='violetred', 'IC4'='purple4', 'IC5'='beige'),
  branch = c('A'=colors_set1[4], 'B'=colors_set1[3], 'C'=colors_set1[2], 'D'=colors_set1[1], 'E'=colors_set1[5])
)

# pseudotime sorted as barplot
ha1 = HeatmapAnnotation(Pseudotime=anno_barplot(mm[sidx,]$pseudotime, height = unit(1, "cm"),
                                                border=F, gp = gpar(fill = "black")),
                        Proliferative_Index=anno_barplot(mm[sidx,]$proliferative_index, height = unit(1, "cm"), 
                                                   border=F, gp = gpar(fill = "black")),
                                                   annotation_name_gp = gpar(fontsize=20));

# features (categorical) - top annotation
features <- c('pre_post', 'icluster', 'subtype_pam50', 'branch_epg20')
anno_df <- data.frame(mm[sidx,features]);
rownames(anno_df) <- mm[sidx,]$sample_id;
colnames(anno_df) <- c('prepost', 'icluster', 'pam50', 'branch');
col_ha = HeatmapAnnotation(df=anno_df, col=annot_colors, 
                           annotation_name_gp = gpar(fontsize=20))

ht1 <- Heatmap(as.matrix(foo1), cluster_columns=F, cluster_rows=T, 
               row_title = "Hallmark",
               row_names_gp = grid::gpar(fontsize=17),
               show_column_names=F, height = nrow(foo1)*unit(5, "mm"), width = ncol(foo1)*unit(2, "mm"),
               name='Hallmark');

ht2 <- Heatmap(as.matrix(foo2), cluster_columns=F, cluster_rows=T, 
               top_annotation = col_ha,
               row_title='Gene',
               row_names_gp = grid::gpar(fontsize=17),
               show_column_names=F, height = nrow(foo2)*unit(5, "mm"), width = ncol(foo2)*unit(2, "mm"),
               name='Gene');

ht3 <- Heatmap(as.matrix(foo3), cluster_columns=F, cluster_rows=T, 
               row_title='Paloma3',
               row_names_gp = grid::gpar(fontsize=17),
               show_column_names=F, height = nrow(foo3)*unit(5, "mm"), width = ncol(foo3)*unit(2, "mm"),
               name='Paloma3');

ht4 <- Heatmap(as.matrix(foo4), cluster_columns=F, cluster_rows=T, 
               row_title='PAM50',
               row_names_gp = grid::gpar(fontsize=17),
               show_column_names=F, height = nrow(foo4)*unit(5, "mm"), width = ncol(foo4)*unit(2, "mm"),
               name='PAM50');

# vertical concatenation of heatmaps
ht_list <- ha1 %v% ht2 %v% ht3 %v% ht1
draw(ht_list, heatmap_legend_side="left", annotation_legend_side="left",
     legend_grouping = "original", legend_labels_gp = gpar(fontsize = 25))

dev.off();






