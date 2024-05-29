#!/usr/bin/env Rscript
#
# Zhengyan 'George' Kan
# May 18, 2023
#
# script to run pre/post comparison for molecular features
#
# tempus_ms_prepost_compare_features.R
#
#
datadir = 'Input';
outdir = 'Output.Manuscript';

library(ggplot2)
library(survival);
library(survminer);
library(lmerTest);
library(EnhancedVolcano)
source('Code/TEMPUS_BC.R');

#feature_type = 'expression';
feature_type = 'molecular';

if (feature_type == 'expression'){

  # gene expression features
  outfile <- paste(outdir, "tempus_ms_prepost_compare_features_continuous.expr.txt", sep='/')
  outfile2 <- paste(outdir, "tempus_ms_prepost_compare_features_discrete.txt", sep='/')  
  plotfile <- paste(outdir, 'tempus_ms_prepost_compare_features.expr.pdf', sep='/');
} else {
  
  # only molecular features (calculated summary)
  outfile <- paste(outdir, "tempus_ms_prepost_compare_features_continuous.txt", sep='/')
  outfile2 <- paste(outdir, "tempus_ms_prepost_compare_features_discrete.txt", sep='/')
  plotfile <- paste(outdir, 'tempus_ms_prepost_compare_features.pdf', sep='/');
}

#----------------------------------------
# (1) load Rdata and molecular data
#
rdata = tempus_bc();
sdata = rdata$sdata;

# sample annotation with covariates
annotation <- sdata[, c('sample_id', 'patient_id', 'prepost', 'tumor_purity', 
                        'organ', 'proliferative_index', 'subcohorts')];
annotation$time <- factor(annotation$prepost, levels=c("Pre", "Post"))

# get molecular features and feature metadata
molfeatfile <- paste(datadir, 'Tempus_molecular_features.txt', sep='/');
mdata <- read.table(molfeatfile, header=T, sep='\t');
metadatafile <- paste(datadir, 'Tempus_molecular_features_metadata.txt', sep='/');
meta <- read.table(metadatafile, header=T, sep='\t');

# align with sample annot
midx <- match(sdata$sample_id, mdata$sample_id)
mdata <- mdata[midx,];

#----------------------------------------
# (2) compare Pre/Post using LMER
# continuous
#

# extract continuous variables
if (feature_type == 'expression'){
  idx <- colnames(mdata) %in% meta$variable[meta$variable_type == 'numerical' &
                                              meta$feature_type %in% c('gene expression')]
} else {
  idx <- colnames(mdata) %in% meta$variable[meta$variable_type == 'numerical' &
                                              meta$feature_type %in% c('molecular')]
}
foo <- mdata[,idx];
rownames(foo) <- mdata$sample_id;

df <- NULL
for(i in 1:ncol(foo)){
  print(colnames(foo)[i])  
  #  fm <- lmer(t(expression[i,])~time+tumor_purity+organ+proliferative_index+(1|patient_id),data=annotation)
  fm <- lmer(foo[,i]~time+tumor_purity+organ+(1|patient_id),data=annotation)
  coef_time <- fixef(fm)[2]
  fval <- anova(fm)["time","F value"]
  pval <- anova(fm)["time","Pr(>F)"]
  means <- tapply(foo[,i], annotation$prepost, mean, na.rm=T);
  #delta <- 100*(means['Post'] - means['Pre'])/mean(foo[,i], na.rm=T);
  delta <- 100*((means['Post'] - means['Pre'])/means['Pre']);
  df = rbind(df, data.frame(variable=colnames(foo)[i], 
                            coef_time=unname(coef_time), 
                            delta=unname(delta),
                            fval=unname(fval), 
                            pval=unname(pval)))
  
}
df$qval <- p.adjust(df$p, method="fdr", n=length(df$p))
df$logp <- ifelse(df$coef_time > 0, -log10(df$pval), log10(df$pval))

write.table(df, file=outfile, quote=F, sep="\t", row.names=F)

#----------------------------------------
# (3) compare Pre/Post using chisq.test
# categorical
#

# extract categorical variables
idx2 <- colnames(mdata) %in% meta$variable[meta$variable_type != 'numerical' &
                                             meta$feature_type %in% c('molecular')]
foo2 <- mdata[,idx2];
rownames(foo2) <- mdata$sample_id;

df2<- NULL
for(i in 1:ncol(foo2)){
  print(colnames(foo2)[i])  
  #res <- wilcox.test(table(annotation$time, foo2[,i]))
  res <- chisq.test(table(annotation$time, foo2[,i]))
  stats <- res$statistic;
  pval <- res$p.value;
  df2 = rbind(df2, data.frame(variable=colnames(foo2)[i], 
                            statistics=unname(stats), 
                            pval=unname(pval)))
}
df2$qval <- p.adjust(df2$p, method="fdr", n=length(df2$p))
write.table(df2, file=outfile2, quote=F, sep="\t", row.names=F)


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
#colors_prepost = palette('Set1')[1:2]
colors_prepost = c('Pre'='lightgrey', 'Post'='grey20');
colors_pam50 = c("purple", "green", "yellow", "orange", "pink");

for (i in 1:length(features)){
#for (i in 1:1){
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
# (5) waterfall plot of hits
#
sidx <- order(df$logp, decreasing=T)
signif <- (df$qval[sidx] < 0.1);
cols <- c('cyan', 'magenta')[as.numeric(signif)+1]
par(mar=c(5,10,4,4))
barplot(df[sidx,]$logp, col=cols, border=F, space=0, ylab='log10(p)', cex.axis=1.5, cex.lab=1.5);
legend("topright", fill=c('magenta'), legend=c('FDR < 10%'), border=F, bty='n', cex=1.5);

#----------------------------------------
# (6) Volcano plot
#
selectedlabs <- c('pam50_cor_her2', 'pam50_cor_luma',
                  'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F10_estrogen_response', 'Paloma3_F8_EMT',
                  'proliferative_index', 'mutsig_13_sigma', 'tmb');
selectedlabs2 <- c('nmf_factor11_proliferation', 'Paloma3_F1_MYC_E2F_activation', 'Paloma3_F9',
                  'pam50_cor_her2', 'proliferative_index');
foo <- df;
EnhancedVolcano(foo, pointSize = 8,
                x='delta', y='pval',
                xlim=c(-200, 300), ylim=c(-1,10),
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
# (7) PAM50 vs. pre/post (Mosaic Plot)
#
par(mfrow=c(1,2))
mdata$prepost <- factor(mdata$prepost, levels=c('Pre', 'Post'))
mosaicplot(table(mdata$prepost, foo2$subtype_pam50), las=1, col=colors_pam50, main='', cex=1.5);

# collect statistics
prop.table(table(mdata$prepost, foo2$subtype_pam50 == 'Her2'), 1)

#--------------------------------------
# (8) PAM50 vs. Pre/Post (Bar plot)
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
