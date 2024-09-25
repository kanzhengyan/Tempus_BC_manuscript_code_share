#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 11, 2024
# Examine tissue site distribution for IC5
#
setwd('../Fig2.FigS7-S8.icluster_characterize')
datadir = 'Input';
outdir = 'Output.Manuscript';

library(RColorBrewer);
library('ggpubr');

plotfile <- paste(outdir, 'FigS7.icluster_tissue_distributions.pdf', sep='/');

#-------------------------------------------------------
# (1) Load updated features, add missing features
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');
mdata$pfs <- mdata$pfs1;
mdata$pfs_event <- mdata$pfs1_censored;
mdata$prepost <- mdata$pre_post;
mdata$tumor_purity <- mdata$tumor_purity_pathologist;

# specifiy color scheme
colors_icluster = c(brewer.pal(4, 'Set1'), 'gray')

#-------------------------------------------------------
# (2) Clean up tissue site annotation
#
tmp = tolower(mdata$tissue_origin);
tmp[grepl('breast', tmp)] = 'breast';
tmp[grepl('liver', tmp)] = 'liver';
tmp[grepl('lymph node', tmp)] = 'lymph node';
tmp[grepl('lung', tmp)] = 'lung';
tmp[grepl('skin', tmp)] = 'skin';
tmp[grepl('bone', tmp) | grepl('vertebral', tmp) | grepl('pelvis', tmp) | grepl('rib', tmp)] = 'bone';
tmp[grepl('pleura', tmp)] = 'pleura';
tmp[grepl('thorax', tmp)] = 'thorax';
tmp[grepl('peritoneum', tmp)] = 'peritoneum';
tmp[grepl('limb', tmp)] = 'limb';
tmp[grepl('abdomen', tmp)] = 'abdomen';
tmp[grepl('cerebellum', tmp) | grepl('lobe', tmp)] = 'peritoneum';
tmp[tmp %in% names(table(tmp))[table(tmp) < 2]] = 'other sites';
mdata$tissue_origin = tmp;

# Graphics output file
pdf(file=plotfile, onefile=T, height=10, width=10, pointsize=12);

#-------------------------------------------------------
# (3) Plot tissue distribution vs. IC1-5
#
y <- 'tissue_origin'
n <- nlevels(as.factor(mdata[,y]));
colors_tissue <- rainbow(n);
#colors_tissue = colorRampPalette(brewer.pal(9, 'Set1'))(n)
foo <- table(mdata[,y]);
foo <- sort(foo, decreasing=T)
par(mar=c(5,6,3,1), mfrow=c(1,2))
barplot(prop.table(table(mdata[,y], mdata[,'icluster']), 2), main='', col=colors_tissue, ylab='% samples', las=2, 
        cex.axis=1.5, cex.names=2, cex.lab=2,
        border='grey');

par(mar=c(5,1,3,3))
frame();
legend('left', legend=levels(as.factor(mdata[,y])), fill=colors_tissue, title='Tissue sites', cex=1.5, bty='n')

#-------------------------------------------------------
# (4) Plot tumor purity vs. IC1-5
#
comps = list(c('IC2', 'IC3'), c('IC3', 'IC4'), c('IC3', 'IC5'));
orders = c('IC1', 'IC2', 'IC3', 'IC4', 'IC5');
ylabel = 'Tumor Purity (%)'
p <- ggboxplot(mdata, 'icluster', 'tumor_purity_tempus', xlab='', ylab=ylabel, title='', outlier.shape=NA, remove=NA,
          scales='free_y', ncol=length(vars),
          fill='icluster', palette=colors_icluster, order=orders) + 
  stat_compare_means(size=5, label='p.signif', comp=comps) + 
  theme(text = element_text(size = 12, face='bold')) + 
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);

ggarrange(p, nrow=2, ncol=3, legend=F)

dev.off()
