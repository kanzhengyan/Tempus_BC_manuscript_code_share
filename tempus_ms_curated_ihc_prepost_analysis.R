#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 12, 2024
# Analyze key IHC marker data from Tempus BC
# Longitudinal changes in Primary vs. Met-Pre vs. Met-Post
#
library(reshape2)
library('ggpubr');
library(RColorBrewer);

datadir = 'Input';
outdir = 'Output.Manuscript';

# input/out files
ihcresfile <- paste(datadir, 'tempus_curated_ihc_results_selected.txt', sep='/');
plotfile <- paste(outdir, 'tempus_ihc_results_analysis.pdf', sep='/');

colors_pam50 = c('Basal'='purple', 'Her2'='green', 'LumA'='yellow', 'LumB'='orange', 'Normal'='pink');
#colors_prepost = c('Pre'='lightgrey', 'Post'='grey12');
#colors_prepost = c('Primary'='white', 'Pre'='lightgrey', 'Post'='grey12');
colors_prepost <- brewer.pal(3, 'Set1')[c(3,2,1)]

# output to text file
res <- read.table(ihcresfile, sep="\t", header=T);

res$Post_subtype <- factor(res$Post_subtype, levels=c('LumA', 'LumB', 'Her2', 'Basal', 'Normal'));
res$Post_icluster <- factor(res$Post_icluster, levels=c('IC1','IC2','IC3','IC4','IC5'));
res$setting <- factor(res$setting, levels=c('Primary', 'Metastatic', 'Other'));

#-------------------------------------------------------------------------
# Compare numerical result vs. Primary/Pre/Post (longitudinal pairs)
#
foo <- res[!is.na(res$numerical_result) & res$prepost_setting %in% c('Primary', 'MetPre', 'MetPost') &
             res$Post_subtype %in% c('LumA', 'LumB', 'Her2'), ];
tmp <- dcast(foo, marker_canonical_name+patient_id~prepost_setting, value.var='numerical_result', mean);

# selected pairs: Pre/Post or Primary/Post
#tmp2<- tmp[!is.na(tmp$MetPre) & !is.na(tmp$MetPost), 1:5];
tmp2<- tmp[(!is.na(tmp$MetPre) & !is.na(tmp$MetPost)) | (!is.na(tmp$Primary) & !is.na(tmp$MetPost)), 1:5];

tmp3 <- melt(tmp2);
colnames(tmp3) <- c('marker', 'patient_id', 'prepost_setting', 'numerical_result') 
dm <- unique(res[,c('patient_id', 'Post_subtype', 'Post_icluster')])
foo_paired <- merge(tmp3, dm, by='patient_id');

foo_paired$prepost_setting <- factor(foo_paired$prepost_setting, levels=c('Primary', 'MetPre', 'MetPost', 'Other'));
foo_paired$Post_subtype <- factor(foo_paired$Post_subtype, levels=c('LumA', 'LumB', 'Her2', 'Basal', 'Normal'));
#foo_paired$Post_icluster2 <- ifelse(foo_paired$Post_icluster == 'IC1', 'IC1', 'IC2-4');
foo_paired$Post_icluster2 <- ifelse(foo_paired$Post_icluster == 'IC1', 'IC1',
                             ifelse(foo_paired$Post_icluster == 'IC4', 'IC4', 'IC2-3'));
markers <- c('ESR1', 'PGR', 'ERBB2')

pdf(file=plotfile, onefile=T, height=8, width=17, pointsize=12);

plist <- list();
c = 1;
for (i in 1:length(markers)){
  marker <- markers[i];
  ymin <- 0;
  ymax <- max(foo_paired[foo_paired$marker == marker,]$numerical_result, na.rm=T)*1.1;
  
  # compare Pre/Post, facet by PAM50 subtype (Post)
  plist[[c]] <- ggpaired(foo_paired[foo_paired$marker == marker,], x='prepost_setting', y='numerical_result', 
                         title=paste(marker, '(Post-PAM50)'), fill='prepost_setting', id='patient_id', ylab='% positive', xlab='', ncol=4,
                         line.color='gray', line.size=.5, point.size = 2,
                         palette=colors_prepost, facet.by='Post_subtype') +
    stat_compare_means(paired=TRUE, size=5, label='p.format', label.y=ymax) +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    theme(text = element_text(size = 14, face='bold'), legend.position = "none");
  c = c + 1;
  
  # compare setting, facet by icluster (Post)
  plist[[c]] <- ggpaired(foo_paired[foo_paired$marker == marker & 
                                      foo_paired$Post_icluster2 %in% c('IC1', 'IC2-3', 'IC4'),], 
#                                      foo_paired$Post_icluster %in% c('IC1', 'IC2', 'IC3', 'IC4'),], 
                         x='prepost_setting', y='numerical_result', 
                         title=paste(marker, '(Post-icluster)'), fill='prepost_setting', id='patient_id', 
                         ylab='% positive', xlab='', ncol=4,
                         line.color='gray', line.size=.5, point.size = 2,
                         palette=colors_prepost, facet.by='Post_icluster2') +
    stat_compare_means(paired=TRUE, size=5, label='p.format', label.y=ymax) +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    theme(text = element_text(size = 14, face='bold'), legend.position = "none");
  c = c + 1;
}

print(ggarrange(plotlist=plist, nrow=2, ncol=2, legend=T, common.legend=T));

dev.off();


