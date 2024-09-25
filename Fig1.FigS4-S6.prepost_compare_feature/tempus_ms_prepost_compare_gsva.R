#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 23, 2023
# Sept. 20, 2024
#
# script to run pre/post comparison for gene signatures
#
setwd('../Fig1.FigS4-S6.prepost_compare_feature')
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2)
library(lmerTest);
library(EnhancedVolcano);

outfile <- paste(outdir, "Table-S4.prepost_compare_gsva.txt", sep='/');
plotfile <- paste(outdir, 'Fig-1d.Fig-S4.prepost_compare_gsva.pdf', sep='/');

#------------------------------
# (1) get sample annotation
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mfdata <- read.table(molfeatfile, header=T, sep='\t');
mfdata$pfs <- mfdata$pfs1;
mfdata$pfs_event <- mfdata$pfs1_censored;
mfdata$prepost <- mfdata$pre_post;
mfdata$tumor_purity <- mfdata$tumor_purity_pathologist;
annotation <- mfdata[, c('sample_id', 'patient_id', 'prepost', 'tumor_purity', 
                         'organ', 'proliferative_index', 'subcohorts')];
annotation$time <- factor(annotation$prepost, levels=c("Pre", "Post"))

# get gene signatures (gsva)

gsvafile <- paste(datadir, 'tempus_ms_supp_table_01.hallmark_gsva.txt', sep='/');
gsvadata <- read.table(gsvafile, sep='\t', header=T)
rownames(gsvadata) <- gsvadata[,1]
gsvadata <- gsvadata[mfdata$sample_id,];

#----------------------------------------
# (2) compare Pre/Post using LMER
# continuous
#
foo <- gsvadata[,-1];
df <- NULL
for(i in 1:ncol(foo)){
  print(colnames(foo)[i])  
  fm <- lmer(foo[,i]~time+tumor_purity+organ+(1|patient_id),data=annotation)
  fmean <- tapply(foo[,i], annotation$prepost, mean);
  delta <- fmean["Post"] - fmean["Pre"];
  coef_time <- fixef(fm)[2]
  fval <- anova(fm)["time","F value"]
  pval <- anova(fm)["time","Pr(>F)"]
  df = rbind(df, data.frame(variable=colnames(foo)[i], 
                            coef_time=unname(coef_time), 
                            fval=unname(fval), 
                            delta=unname(delta),
                            pval=unname(pval)))
  
}
df$qval <- p.adjust(df$p, method="fdr", n=length(df$p))
df$logp <- ifelse(df$coef_time > 0, -log10(df$pval), log10(df$pval))
write.table(df, file=outfile, quote=F, sep="\t", row.names=F)

#----------------------------------------
# (3) boxplot Pre/Post comparison
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=10, pointsize=12);

pcut <- 0.05;
df <- df[order(df$pval),];
features <- df[df$pval < pcut,]$variable;

# color labels
colors_organ = rainbow(4);
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');
colors_pam50 = c("purple", "green", "yellow", "orange", "pink");

features <- c('E2F_TARGETS', 'MYC_TARGETS_V1')
for (i in 1:length(features)){
  feature <- features[i];
  x <- foo[,feature];

  c = 1;
  plist <- list();
    
  # (1) Examine variable change vs. Pre/Post
  #
  comps <- list(c("Pre", "Post"));
  mdata = data.frame(feature=x, annotation);
  ymax <- max(mdata$feature);
  ymin <- min(mdata$feature);
  plist[[c]] = ggboxplot(mdata, "prepost", "feature", xlab='', ylab='GSVA', 
                         title='', outlier.shape=NA,
                         ylim=c(ymin, ymax),
                         fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post')) +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    stat_compare_means(comparison=comps, size=4, label='p.format') +
    theme(text = element_text(size = 12, face='bold'), legend.position='none') +
    geom_jitter(pch=21, fill='gray', width=0.2, alpha=.6, size=3);
  c = c+1;
  
  # (2) Examine longitudial change vs. Pre/Post
  #
  mdata_paired = mdata[mdata$subcohorts == '1-Paired',]
  plist[[c]] <- ggpaired(mdata_paired, x='prepost', y='feature', title='', fill='prepost', 
                         id='patient_id', xlab='', ylab='',
                         ylim=c(ymin, ymax),
                         line.color='gray', line.size=.4, palette=colors_prepost, order=c('Pre', 'Post')) +
    scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
    stat_compare_means(paired=TRUE, size=4, label='p.format') +
    theme(text = element_text(size = 12, face='bold'), legend.position = "none");
  c = c+1;

  p <- ggarrange(plotlist=plist, nrow=2, ncol=4);
  print(annotate_figure(p, fig.lab=feature, fig.lab.size=12, fig.lab.face='bold'))
}

#----------------------------------------
# (4) Volcano plot
#
selectedlabs <- c('GLYCOLYSIS', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'PEROXISOME',
                  'MTORC1_SIGNALING', 'OXIDATIVE_PHOSPHORYLATION');
selectedlabs2 <- c('G2M_CHECKPOINT', 'GLYCOLYSIS', 'MTORC1_SIGNALING', 'MYC_TARGETS_V1', 
                   'PI3K_AKT_MTOR_SIGNALING', 'UNFOLDED_PROTEIN_RESPONSE', 'UV_RESPONSE_UP');

foo <- df;
EnhancedVolcano(foo, pointSize = 8,
                x='coef_time', y='pval',
                xlim=c(-.3, .3), ylim=c(0,13),
                FCcutoff=.1, pCutoff=0.05,
                lab=NA,
#                lab=foo$variable,
#                selectLab = c('UNFOLDED_PROTEIN_RESPONSE'),
                boxedLabels = T, labFace='bold', 
                drawConnectors=TRUE, colConnectors='black', widthConnectors=1,
                xlab = 'Coef',
                gridlines.major = F, gridlines.minor = F,
                legendLabels = c("NS", "Coef", "p-value", expression("p-value" ~ and ~ "Coef")),)

dev.off();
