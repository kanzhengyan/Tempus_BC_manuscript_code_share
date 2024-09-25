#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# August 16, 2023
# Sept. 20, 2023
#
# script to characterize gene expression and gene signatures
# using multi-panel boxplots
#
setwd('../Fig3.FigS9.cse_gsva_DE_analysis')
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(RColorBrewer);

plotfile <- paste(outdir, 'Fig-3.cse_gsva_multipanel.pdf', sep='/');
plotfile2 <- paste(outdir, 'Fig-S9.cse_gene_multipanel.pdf', sep='/');

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
csefile <- paste(datadir, 'tempus_ms_supp_data.gene_cse.txt', sep='/');
csedata <- read.table(csefile, sep='\t', header=T)
rownames(csedata) <- csedata[,1]
cse <- t(csedata[,-1]);

# load CSE GSVA data
gsvafile <- paste(datadir, 'tempus_ms_supp_data.hallmark_cse_gsva.txt', sep='/');
gsvadata <- read.table(gsvafile, sep='\t', header=T)
rownames(gsvadata) <- gsvadata[,1]
gsva <- t(gsvadata[,-1]);

# add gene signature GSVA data
df2 <- data.frame(sample_id=colnames(gsva), t(gsva));
mdata <- merge(mdata, df2, by='sample_id', all.x=T, sort=F);
colnames(mdata) <- gsub('c2_cgp.', '', colnames(mdata))
colnames(mdata) <- gsub('EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT', colnames(mdata))

#=============================================
# (3) prepare data for analysis
#
pdf(file=plotfile, onefile=T, height=8, width=11, pointsize=12);

# select data variables for analysis  
gsvas_ic1 <- c('ESTROGEN_RESPONSE_EARLY', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'GLYCOLYSIS');
gsvas_ic2 <- c('MTORC1_SIGNALING', 'IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE', 'EMT');
gsvas_her2 = c('SMID_BREAST_CANCER_ERBB2_UP', 'GLYCOLYSIS', 'MTORC1_SIGNALING');
gsvas_basal = c('IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE', 'EMT');
list_markers <- list(gsvas_ic1, gsvas_ic2);

# color labels
annot_colors = list(
  pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink'),
  icluster =  c(brewer.pal(4, 'Set1'), 'grey'),
  prepost = c('Post'='grey20', 'Pre'='lightgrey'),
  mut = c('0'='deepskyblue2', '1'='firebrick2'),
  amp = c('0'='deepskyblue2', '1'='orange3'),
  del = c('0'='deepskyblue2', '1'='darkblue'),
  histology = rainbow(3),
  organ = rainbow(4)
)

#---------------------------------------------------------
# (3.1) compare variable vs. key categories 
# Pre/Post, Pre/During/Post, Tissue, PAM50, IC, 
#
ylabel <- 'GSVA';
for (k in 1:length(list_markers)){
  markers <- list_markers[[k]];
  plist <- list();
  categories <- c('subtype_pam50', 'icluster');
  for (i in 1:length(categories)){
    category <- categories[i];
    
    # append data for multiple variables
    for (j in 1:length(markers)){
      marker <- markers[j];
      df <- data.frame(mdata, value=mdata[,marker], varlabel=marker);
      if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
    }
    
    # reorder levels
    foo$varlabel = factor(foo$varlabel, levels=markers)  
    
    if (category == 'subtype_pam50') {
      comps <- list(c('Basal', 'Her2'), c('Her2', 'LumA'), c('Her2', 'LumB'));
      orders <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
      colors <- annot_colors$pam50;
      foo <- foo[foo$subtype != 'Normal',];
    } else if (category == 'icluster') {
      comps <- list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'), c('IC1', 'IC5'));
      orders <- c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
      colors <- annot_colors$icluster;
      foo <- foo[foo$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4', 'IC5'),]
    }
    nrows = ifelse(length(markers) > 4, 2, 1)
    
    # Luminal markers
    plist[[i]] = ggboxplot(foo, category, 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA,
                           facet.by='varlabel', scales='free_y', nrow=nrows, mult=.5,
                           fill=category, palette=colors,  order=orders) +
      stat_compare_means(comparison=comps, size=5, label='p.signif') +
      scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
      theme(text = element_text(size = 12, face='bold')) +
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2)
    
  }
  print(ggarrange(plotlist=plist, nrow=2, ncol=1, legend=F));
}

#---------------------------------------------------------
# (3.2) compare variable vs. Pre/Post (paired)
#
for (k in 1:length(list_markers)){
  markers <- list_markers[[k]];

  # data for only sample pairs
  mdata_paired = mdata[mdata$subcohorts == '1-Paired',]
  
  # append data for multiple variables
  for (j in 1:length(markers)){
    marker <- markers[j];
    df <- data.frame(mdata_paired, value=mdata_paired[,marker], varlabel=marker);
    if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
  }
  foo$varlabel = factor(foo$varlabel, levels=markers);
  
  p2 <- ggpaired(foo, x='prepost', y='value', title='', fill='prepost', id='patient_id', xlab='', ylab=ylabel,
                facet.by='varlabel', scale='free_y', ncol=length(markers),
                line.color='gray', line.size=.4, palette=annot_colors$prepost, order=c('Pre', 'Post')) +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    stat_compare_means(paired=TRUE, size=4, label='p.format') +
    theme(text = element_text(size = 12, face='bold'), legend.position = "none");

  print(ggarrange(p2, nrow=2, ncol=1, widths=c(4,1)));
}

dev.off();

#=====================================================================

# add gene CSE data
df1 <- data.frame(sample_id=colnames(cse), t(cse));
mdata <- merge(mdata, df1, by='sample_id', all.x=T, sort=F);

pdf(file=plotfile2, onefile=T, height=8, width=11, pointsize=12);

# select data variables for analysis
markers_luminal <- c('FOXA1', 'SPDEF', 'KRT18');
markers_basal <- c('ERBB2', 'KRT17', 'KRT14', 'KRT5');
markers_emt <- c('BMP1', 'SNAI2', 'VIM');
list_markers <- list(markers_luminal, markers_basal, markers_emt);

ylabel <- 'Log2(TPM)';
for (k in 1:length(list_markers)){
  markers <- list_markers[[k]];
  plist <- list();
  categories <- c('subtype_pam50', 'icluster');
  for (i in 1:length(categories)){
    category <- categories[i];
    
    # append data for multiple variables
    for (j in 1:length(markers)){
      marker <- markers[j];
      df <- data.frame(mdata, value=mdata[,marker], varlabel=marker);
      if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
    }
    
    # reorder levels
    foo$varlabel = factor(foo$varlabel, levels=markers)  
    
    if (category == 'subtype_pam50') {
      comps <- list(c('Basal', 'Her2'), c('Her2', 'LumA'), c('Her2', 'LumB'));
      orders <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
      colors <- annot_colors$pam50;
      foo <- foo[foo$subtype != 'Normal',];
    } else if (category == 'icluster') {
      comps <- list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'), c('IC1', 'IC5'));
      orders <- c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
      colors <- annot_colors$icluster;
      foo <- foo[foo$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4', 'IC5'),]
    }
    nrows = ifelse(length(markers) > 4, 2, 1)
    
    # Luminal markers
    plist[[i]] = ggboxplot(foo, category, 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA,
                           facet.by='varlabel', scales='free_y', nrow=nrows, mult=.5,
                           fill=category, palette=colors,  order=orders) +
      stat_compare_means(comparison=comps, size=5, label='p.signif') +
      scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
      theme(text = element_text(size = 12, face='bold')) +
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2)
    
  }
  print(ggarrange(plotlist=plist, nrow=2, ncol=1, legend=F));
}

dev.off();
