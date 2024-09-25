#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 18, 2023
# Sept. 20, 2024
#
# script to run pre/post comparison for molecular features
#
setwd('../Fig1.FigS4-S6.prepost_compare_feature')
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2);
library(ggpubr);
library(lmerTest);
library(EnhancedVolcano);

outfile <- paste(outdir, "Table-S2.prepost_compare_features.txt", sep='/')
plotfile <- paste(outdir, 'Fig-1c.Fig-S4-S6.prepost_compare_features.pdf', sep='/');

#---------------------------------------------------------
# (1) get molecular features and feature metadata
#
molfeatfile <- paste(datadir, 'Tempus_molecular_features_extended.txt', sep='/');
mfdata <- read.table(molfeatfile, header=T, sep='\t');
metadatafile <- paste(datadir, 'Tempus_molecular_features_extended_metadata.txt', sep='/');
meta <- read.table(metadatafile, header=T, sep='\t');

mfdata$pfs <- mfdata$pfs1;
mfdata$pfs_event <- mfdata$pfs1_censored;
mfdata$prepost <- mfdata$pre_post;
mfdata$tumor_purity <- mfdata$tumor_purity_pathologist;
annotation <- mfdata[, c('sample_id', 'patient_id', 'prepost', 'tumor_purity', 
                        'organ', 'proliferative_index', 'subcohorts')];
annotation$time <- factor(annotation$prepost, levels=c("Pre", "Post"))

#----------------------------------------
# (2) compare Pre/Post using LMER
# continuous variables
#

# extract continuous variables
idx <- colnames(mfdata) %in% meta$variable[meta$variable_type == 'numerical' & 
                                             meta$feature_type %in% c('molecular')]

foo <- mfdata[,idx];
rownames(foo) <- mfdata$sample_id;

df1 <- NULL
for(i in 1:ncol(foo)){
  print(colnames(foo)[i])  
  fm <- lmer(foo[,i]~time+tumor_purity+organ+(1|patient_id),data=annotation)
  coef_time <- fixef(fm)[2]
  fval <- anova(fm)["time","F value"]
  pval <- anova(fm)["time","Pr(>F)"]
  means <- tapply(foo[,i], annotation$prepost, mean, na.rm=T);
  delta <- 100*((means['Post'] - means['Pre'])/abs(means['Pre']));
  df1 = rbind(df1, data.frame(variable=colnames(foo)[i], 
                            coef_time=unname(coef_time), 
                            delta=unname(delta),
                            fval=unname(fval), 
                            pval=unname(pval)))
  
}
df1$qval <- p.adjust(df1$p, method="fdr", n=length(df1$p))
df1$logp <- ifelse(df1$coef_time > 0, -log10(df1$pval), log10(df1$pval))

# extract categorical variables
idx2 <- colnames(mfdata) %in% meta$variable[meta$variable_type != 'numerical' &
                                             meta$feature_type %in% c('molecular')]
foo2 <- mfdata[,idx2];
rownames(foo2) <- mdata$sample_id;

foo2[foo2 == 'Yes'] <- 1;
foo2[foo2 == 'No'] <- 0;
foo2$brca_pathogenic_mutation <- as.numeric(foo2$brca_pathogenic_mutation)
foo2$bc_susceptibility_gene_pathogenic_mutation <- as.numeric(foo2$bc_susceptibility_gene_pathogenic_mutation)

df2<- NULL
for(i in 1:ncol(foo2)){
  print(colnames(foo2)[i])
  if (colnames(foo2)[i] %in% c('subtype_pam50', 'mutsig_sigma_status')){
    res <- chisq.test(table(annotation$time, foo2[,i]))
    stats <- res$statistic;
    lfc <- 0;
  } else {
    res <- fisher.test(table(annotation$time, foo2[,i]));
    stats <- res$estimate;
    tmp <- prop.table(table(annotation$time, foo2[,i]), 1);
    lfc <- log2(tmp['Post','1']/tmp['Pre','1']);
  }
  pval <- res$p.value;
  df2 = rbind(df2, data.frame(variable=colnames(foo2)[i], 
                              statistics=unname(stats), 
                              lfc=lfc,
                              pval=unname(pval)))
}
df2$qval <- p.adjust(df2$p, method="fdr", n=length(df2$p))
df2$logp <- ifelse(df2$lfc > 0, -log10(df2$pval), log10(df2$pval))

#--------------------------------------------------------------------
# write to combined file, both categorical & numerical variables
#
df <- rbind(data.frame(variable=df1$variable, type='numerical', df1[,c('coef_time', 'delta', 'pval', 'logp')]),
            data.frame(variable=df2$variable, type='categorical', coef_time=NA, delta=NA, pval=df2$pval, logp=df2$logp))
df$qval <- p.adjust(df$p, method="fdr", n=length(df$p))

write.table(df, file=outfile, quote=F, sep="\t", row.names=F)

#----------------------------------------
# (4) boxplot Pre/Post comparison
#

# graphics output file
pdf(file=plotfile, onefile=T, height=9, width=10, pointsize=12);

pcut <- 0.05;
df <- df[order(df$pval),]
features <- df[df$pval < pcut,]$variable;

# color labels
colors_organ = rainbow(4);
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');
colors_pam50 = c("purple", "green", "yellow", "orange", "pink");
features <- c('proliferative_index', 'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response')
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
  plist[[c]] = ggboxplot(mdata, "prepost", "feature", xlab='', ylab='', 
                         title='', outlier.shape=NA,
                 fill='prepost', palette=colors_prepost,  order=c('Pre', 'Post'),
                 ylim=c(ymin, ymax)) +
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
# (5) Volcano plot
#
selectedlabs <- c('pam50_cor_her2', 'pam50_cor_luma',
                  'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT',
                  'proliferative_index', 'mutsig_13_sigma', 'tmb');
selectedlabs2 <- c('nmf_factor11_proliferation', 'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F9',
                  'pam50_cor_her2', 'proliferative_index');
foo <- df1;
EnhancedVolcano(foo, pointSize = 8,
                x='delta', y='pval',
                xlim=c(-200, 200), ylim=c(-1,10),
                FCcutoff=20, pCutoff=0.05,
                lab=NA,
                #lab=foo$variable,
                #selectLab = selectedlabs2,
                #drawConnectors=TRUE, colConnectors='black', widthConnectors=1,
                #boxedLabels = T, labFace='bold', 
                xlab = '% Change',
                gridlines.major = F, gridlines.minor = F,
                legendLabels = c("NS", "% Change", "p-value", expression("p-value" ~ and ~ "% Change")),)

#----------------------------------------
# (6) PAM50 vs. pre/post (Mosaic Plot)
#
par(mfrow=c(1,2))
mdata$prepost <- factor(mdata$prepost, levels=c('Pre', 'Post'))
mosaicplot(table(mdata$prepost, foo2$subtype_pam50), las=1, col=colors_pam50, main='', cex=1.5);

#--------------------------------------
# (7) PAM50 vs. Pre/Post (Bar plot)
#
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');

tbl <- prop.table(table(foo2$subtype_pam50, mdata$prepost), 2);
df <- as.data.frame(tbl);
colnames(df) <- c('pam50', 'prepost', 'proportion');
df$prepost <- factor(df$prepost, levels=c('Pre', 'Post'));
df <- df[df$pam50 != 'Normal',]

p <- ggplot(df, aes(x = pam50, y = proportion, fill = prepost)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=colors_prepost) +
  ylab('% Samples') + scale_y_continuous(labels = scales::percent) +
  theme(text = element_text(size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank());

print(ggarrange(p, nrow=2, ncol=2))

dev.off();



