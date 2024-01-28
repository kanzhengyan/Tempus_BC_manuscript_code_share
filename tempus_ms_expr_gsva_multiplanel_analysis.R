#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# August 16, 2023
#
# script to characterize gene expression and gene signatures
# using multi-panel boxplots 
#
# tempus_ms_expr_gsva_multipanel_analysis.R > tempus_ms_expr_gsva_multipanel_analysis.pdf
#
# - extract a list of gene expression (CSE) and GSVA signatures (CSE), plot
#   (1) score vs. Pre/Post time (boxplot, p-value)
#   (2) score vs. Pre/Post time (paired boxplot)
#   (3) score vs. Key categories (IC, IC-PAM50)
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library('ggpubr');
library(RColorBrewer);
source('Code/TEMPUS_BC.R');

plotfile <- paste(outdir, 'tempus_ms_expr_gsva_multipanel_analysis.pdf', sep='/');

#--------------------------------
# Tempus BC R data object
#
rdata = tempus_bc();
sdata = rdata$sdata;

#---------------------------------------------
# (1) get sample metadata
#
mfdatafile <- paste(datadir, "Tempus_molecular_features_extended_v3.txt", sep='/');
mfdata <- read.table(mfdatafile, sep="\t", header=T);

# merge sample annotation w/ molecular features
mdata <- merge(mfdata, sdata[,c('sample_id', 'subcohorts', 'time', 'organ', 'subcohorts')],
               by='sample_id', all.x=T, sort=F);
mdata$prepost = factor(mdata$pre_post, levels=c('Pre', 'Post'))
mdata <- subset(mdata, select=c('sample_id', 'patient_id', 'subcohorts', 'prepost', 'time', 'subtype_pam50', 'icluster'));

#-------------------------------------------
# add icluster+pam50 subgroups
#
groups <- paste(mdata$icluster, mdata$subtype_pam50, sep=':')
selected_groups <- names(table(groups))[table(groups) >= 7];
idx2 <- groups %in% selected_groups & mdata$icluster != 'IC5' & !is.na(mdata$icluster)
mdata$ic_pam50 <- groups;
mdata$ic_pam50[!idx2] <- NA;

#------------------------------------------------------
# (2) add CSE expression & gsva data
#

# cancer-specific expression
csedata = rdata$Expression$BayesPrism_cancer

# CSE gene signatures (GSVA)
gsvadata_cse = rdata$Analysis$GSVA_BayesPrism_cancer;

# add gene CSE data
genes <- c('FOXA1', 'SPDEF', 'KRT18', 'KRT8', 'KRT17', 'KRT14', 'KRT5', 'PRF1', 
                'IFNAR2', 'IL6', 'FAS', 'CD40', 'BMP1', 'SNAI2', 'VIM');
exp <- csedata[genes,];
df1 <- data.frame(sample_id=colnames(exp), t(exp));
mdata <- merge(mdata, df1, by='sample_id', all.x=T, sort=F);

# add gene signature GSVA data
gsvas <- c('h.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'h.HALLMARK_IL6_JAK_STAT3_SIGNALING',
           'h.HALLMARK_INFLAMMATORY_RESPONSE');
gsva <- gsvadata_cse[gsvas,]
df2 <- data.frame(sample_id=colnames(gsva), t(gsva));
mdata <- merge(mdata, df2, by='sample_id', all.x=T, sort=F);
colnames(mdata) <- gsub('h.HALLMARK_', '', colnames(mdata))
colnames(mdata) <- gsub('EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT', colnames(mdata))

#=============================================
# (3) prepare data for analysis
#
pdf(file=plotfile, onefile=T, height=8, width=10, pointsize=12);

# select data variables for analysis  
markers_luminal <- c('FOXA1', 'SPDEF', 'KRT18');
markers_basal <- c('KRT17', 'KRT14', 'KRT5');
markers_immune <- c('PRF1', 'IFNAR2', 'IL6', 'FAS');
markers_emt <- c('BMP1', 'SNAI2', 'VIM');
gsvas_basal = c('IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE', 'EMT');
list_markers <- list(markers_luminal, markers_basal, markers_immune, markers_emt, gsvas_basal);

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

for (k in 1:length(list_markers)){

  markers <- list_markers[[k]];
  if (k == 5) { ylabel = 'GSVA'} else { ylabel = 'Log2(TPM)'}
  
  #---------------------------------------------------------
  # (3.1) compare variable vs. key categories 
  # Pre/Post, Pre/During/Post, Tissue, PAM50, IC, 
  # MUT, AMP
  #
  plist <- list();
  categories <- c('prepost', 'subtype_pam50', 'icluster', 'time');
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
    
    if (category == 'prepost') {
      comps <- list(c("Pre", "Post"));
      orders <- c('Pre', 'Post');
      colors <- annot_colors$prepost;
    } else if (category == 'time'){
      comps <- list(c("BL", "OT"), c("OT", "PD"), c("BL", "PD"));
      orders <- c('BL', 'OT', 'PD');
      colors <- annot_colors$time;
    } else if (category == 'subtype_pam50') {
      comps <- list(c('Basal', 'Her2'), c('Basal', 'LumA'), c('Basal', 'LumB'));
      orders <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal');
      colors <- annot_colors$pam50;
      foo <- foo[foo$subtype != 'Normal',];
    } else if (category == 'icluster') {
      comps <- list(c('IC1', 'IC2'), c('IC1', 'IC3'), c('IC1', 'IC4'));
      orders <- c('IC1', 'IC2', 'IC3', 'IC4');
      colors <- annot_colors$icluster;
      foo <- foo[foo$icluster %in% c('IC1', 'IC2', 'IC3', 'IC4'),]
    }
    
    # Luminal markers
    plist[[i]] = ggboxplot(foo, category, 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA,
                           facet.by='varlabel', scales='free_y', ncol=length(markers), mult=.5,
                           fill=category, palette=colors,  order=orders) +
      stat_compare_means(comparison=comps, size=5, label='p.signif') +
      scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
      theme(text = element_text(size = 12, face='bold')) +
#      theme(text = element_text(size = 9, face='bold')) +
      geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2)
    
  }
  print(ggarrange(plotlist=plist, nrow=2, ncol=1));
    
  #-------------------- 
  # (3.2) IC-PAM50 
  #
  
  # append data for multiple variables
  for (j in 1:length(markers)){
    marker <- markers[j];
    df <- data.frame(mdata, value=mdata[,marker], varlabel=marker);
    if (j == 1) { foo <- df; } else { foo <- rbind(foo, df)}
  }
  foo$varlabel = factor(foo$varlabel, levels=markers)  
  
  comps <- list(c('IC1:Basal', 'IC1:Her2'), c('IC1:Her2', 'IC1:LumB'), c('IC1:Basal', 'IC1:LumB'),
                c('IC1:LumB', 'IC2:LumB'), c('IC1:LumB', 'IC3:LumB'));
  orders <- c('IC1:Basal',  'IC1:Her2', 'IC1:LumB', 'IC2:LumB', 'IC3:LumB');
  colors <- annot_colors$ic_pam50;
  foo <- foo[foo$ic_pam50 %in% c('IC1:Basal', 'IC1:Her2', 'IC1:LumB', 'IC2:LumB', 'IC3:LumB'),]
  
  category = 'ic_pam50';
  p1 = ggboxplot(foo, category, 'value', xlab='', ylab=ylabel, title='', outlier.shape=NA,
                facet.by='varlabel', scales='free_y', ncol=length(markers), mult=.5,
                fill=category, palette=colors,  order=orders) +
    stat_compare_means(comparison=comps, size=5, label='p.signif') +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    theme(text = element_text(size = 12, face='bold')) +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=2) +
    rotate_x_text(45)
  
  print(ggarrange(p1, nrow=2, ncol=1, heights=c(2,1)));

  #---------------------------------------------------------------------------
  # (4) Longitudinal variable changes, Pre vs. Post (paired BL/PT)
  #

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
    stat_compare_means(paired=TRUE, size=5, label='p.format') +
    #theme(text = element_text(size = 12, face='bold'), legend.position = "none");
    theme(text = element_text(size = 9, face='bold'), legend.position = "none");

  #print(ggarrange(p2, nrow=2, ncol=1));  
  print(ggarrange(p2, nrow=2, ncol=2, widths=c(3,2)));
}

dev.off();
